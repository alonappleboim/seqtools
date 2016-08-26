import os
import sys
import numpy
from scipy import io as sio
from abc import ABCMeta, abstractmethod


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
    def export(self, sample_stats):
        '''
        export stats and counts to files in out_path
        :param sample_stats: A dictionary of statistic names, and sample-based dictionary of statistic dictionaries.
        :return: generated file names
        '''
        pass


class TabExporter(Exporter):
    name = 'tab'
    description = 'export to tab delimited files'
    args = {}

    def export(self, sample_features, sample_stats):
        fnames = []
        for name, sdict, stats in sample_stats:
            fname = name + '.tab'
            fnames.append(fname)
            fout = open(self.out_path + os.sep + fname, 'w')
            # header:
            for f in sample_features.values():
                fout.write(str(f) + '\t' + '\t'.join(str(s.fvals[f]) for s in sdict.keys()) + '\n')
            # counts
            for stat in stats:
                fout.write(stat + '\t' + '\t'.join(str(sdict[s][stat]) for s in sdict.keys()) + '\n')
            fout.close()
        return fnames


class MatExporter(Exporter):
    name = 'mat'
    description = 'export to a single matlab file'
    args = {'r': (bool, True, 'whether to include a multidimensional array version of the data')}

    def export(self, stats, counts):
        return []


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