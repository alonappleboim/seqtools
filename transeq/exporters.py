import os
import sys
import numpy as np
from scipy import io as sio
from abc import ABCMeta, abstractmethod
from collections import OrderedDict

def collect_exporters():
    exporters = {}
    exporter_module = sys.modules[__name__]
    for x in dir(exporter_module):
        try:
            if issubclass(exporter_module.__dict__[x], Exporter):
                ename = exporter_module.__dict__[x].name
                if ename is not None:
                    exporters[ename] = exporter_module.__dict__[x]
        except TypeError: pass
    return exporters


class Exporter(object):
    """
    writes outputs to a specific format
    """
    __metaclass__ = ABCMeta
    name = None
    description = None
    args = None

    def __init__(self, out_path, **kwargs):
        self.out_path = out_path
        self.__dict__.update(kwargs)

    @abstractmethod
    def export(self, features, sample_stats):
        """
        export statistics to output_dir

        :param features: the feature collection associated with the statistics
        :param samples: a list of samples, determines output order of all stats
        :param stats: a list of 3-tuples: name - used for identyifing the exported data, sdict - the statistics
        dictionary, indexed by sample, and then by statistic, and stats - the statistic ordered names.
        :return: generated file names
        """
        pass


class TabExporter(Exporter):
    name = 'tab'
    description = 'export to tab delimited files'
    args = {}

    def export(self, features, samples, sample_stats):
        fnames = []
        for name, sdict, stats in sample_stats:
            fname = name + '.tab'
            fnames.append(fname)
            fout = open(self.out_path + os.sep + fname, 'w')
            # header:
            for f in features:
                fout.write(str(f) + '\t' + '\t'.join(str(s.fvals[f]) for s in samples) + '\n')
            # counts
            for stat in stats:
                fout.write(stat + '\t' + '\t'.join(str(sdict[s][stat]) for s in samples) + '\n')
            fout.close()
        return fnames


class MatExporter(Exporter):
    name = 'mat'
    description = 'export to a single matlab file'
    args = {'r': (bool, True, 'whether to include a multidimensional array version of the data'),
            'name': (str, 'tts', 'used as the file name and the matlab struct name.')}

    def export(self, features, samples, sample_stats):
        fpath = self.out_path + os.sep + self.name
        lg = {}
        for f in features:
            dtype = np.object if f.type is str else np.double
            lg[f.name] = np.array([s.fvals[f] for s in samples], dtype=dtype)

        s = dict(lg=lg)
        for name, sdict, stats in sample_stats:
            lg[name] = np.array(stats, dtype=np.object)
            s[name] = np.array([[float(sdict[s][stat]) for s in samples] for stat in stats])

        if self.r:
            lg = OrderedDict()
            r = {}
            for f in features: lg[f.name] = sorted(f.vals)
            for name, _, _ in sample_stats:
                dims = (s[name].shape[0],) + tuple([len(f.vals) for f in features])
                arr = np.empty(dims)
                for si, arri in enumerate(MatExporter.sample2arr(samples, features, lg)):
                    itup = (slice(None),) + arri
                    arr[itup] = s[name][:,si]
                r[name] = arr
            for f in features:
                dtype = np.object if f.type is str else np.double
                lg[f.name] = np.array(sorted(f.vals), dtype=dtype)
            r['lg'] = lg
            s['r'] = r
        s['doc'] = np.array('"lg" is a short for "legend", there is a complex one ("lg") for the samples'
                            ' and then another one for the rows of each table. "r_" is a prefix for reshaped '
                            'data or respective legends',dtype=str)
        sio.savemat(fpath, {self.name: s})
        return [self.name+'.mat']

    @staticmethod
    def sample2arr(samples, features, rlg):
        for s in samples:
            r_idx = tuple(rlg[fname].index(s.fvals[f]) for f, fname in zip(features, rlg.keys()))
            yield r_idx


class NumpyExporter(Exporter):
    name = 'np'
    description = 'pickle numpy versions of the data'
    args = {'r': (bool, True, 'whether to include a multidimensional array version of the data')}

    def export(self, stats, counts):
        pass


def exporters_from_string(estring, out_path):
    eclasses = collect_exporters()
    es = []
    for e, argdict in parse_exporters(estring):
        if e not in eclasses:
            raise ValueError('No exporter %s' % e)
        E = eclasses[e]
        args = {}
        for aname, aval in argdict.items():
            args[aname] = E.args[aname][0](aval)
        for a, props in E.args.items():
            if a not in args: args[a] = props[1]  # default
        es.append(E(out_path, **args))
    return es


def parse_exporters(estr):
    exps = []
    for es in estr.split(';'):
        name, arg_s = es.split('(')
        arg_s = arg_s.strip()[:-1].split(',')
        args = {}
        for a in arg_s:
            a = a.strip()
            if not a: continue
            aname, aval = a.split('=')
            args[aname] = aval
        exps.append((name,args))
    return exps