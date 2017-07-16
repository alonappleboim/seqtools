#! /cs/bd/tools/nflab_env/bin/python3.4
import sys
import os
import re
import scipy.io as sio
import numpy as np
from scipy import sparse
import argparse
from collections import OrderedDict

# INTERPRETER = '/cs/bd/tools/nflab_env/bin/python3.4'
# if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
#     import subprocess as sp
#     import os
#     scriptpath = os.path.abspath(sys.modules[__name__].__file__)
#     sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
#     exit()


def get_data(path, var_parsers):
    """
    Parse the files in path in light of given regexps

    :param path: a path containing files corresponding to given regexps
    :param var_parsers: a list of regexps that capture variables
    :return:
        vars - an ordered dictionary of sample variables, mapping to ordered lists of variable values
        file_map - a map from variable value tuples to matching .mat file paths
        file_list - a list of files that matched criteria, in order of iteration
        var_lists - a list of lists, one per variable, of variable values, corresponding to files in file_list
    """
    file_map, vars = {}, OrderedDict()
    file_list, var_lists = [], []
    for vp in var_parsers:
        for v in vp.groupindex.keys():
            if v not in vars:
                vars[v] = []
                var_lists.append([])
    for f in os.listdir(os.path.abspath(path)):
        for vp in var_parsers:
            m = vp.match(f)
            if m is not None: break
        if m is None: continue  # no match
        fpath = path + os.path.sep + f
        s_dict, s_vals = m.groupdict(), []
        for i, v in enumerate(vars):
            val = s_dict[v] if v in s_dict else None
            var_lists[i].append(val)
            s_vals.append(val)
        k = tuple(s_vals)
        if k in file_map:
            raise ValueError('could not differentiate between two samples: %s and %s' % (file_map[k], fpath))
        for var, val in s_dict.items():
            if val not in vars[var]: vars[var].append(val)
        file_list.append(fpath)
        file_map[k] = fpath
    return vars, file_map, file_list, var_lists


def build_matlab_objects(vars, file_map, file_list, var_lists):
    """

    :param vars:
    :param file_map:
    :param file_list:
    :param var_lists:
    :return:
        mdict - a Matlab struct with following fields:
                * d - a cell array of #chromosomes, each holding a matrix of N_ixS, where N_i is the chromosome
                      length, and S is the number of samples.
                * l - a legend for the data, a struct, with following fields:
                      * samples - a list of sample names in the order they appear in "d{i}"
                      * vars - a struct with every field corresponding to a variable of the samples, mapping to a list
                               of values of that variable in the order of samples in "d"
                      * r - a struct with variables reshaped to a multi-dimensional array for easy access to subsets of
                            samples.
                * r - a mapping of indices from the N_ixS structure of the data to multidimensional array, with every
                      dimension corresponding to a different variable, for for easy access to subsets of samples.
    """
    for
    if self.r:
        lg = OrderedDict()
        r = {}
        for f in features: lg[f.name] = sorted(f.vals)
        for name, _, _ in sample_stats:
            dims = (s[name].shape[0],) + tuple([len(f.vals) for f in features])
            arr = np.empty(dims)
            for si, arri in enumerate(MatExporter.sample2arr(samples, features, lg)):
                itup = (slice(None),) + arri
                arr[itup] = s[name][:, si]
            r[name] = arr
        for f in features:
            dtype = np.object if f.type is str else np.double
            lg[f.name] = np.array(sorted(f.vals), dtype=dtype)
        r['lg'] = lg
        s['r'] = r
    s['doc'] = np.array('"lg" is a short for "legend", there is a complex one ("lg") for the samples'
                        ' and then another one for the rows of each table. "r_" is a prefix for reshaped '
                        'data or respective legends', dtype=str)
    sio.savemat(fpath, {self.name: s})
    return [self.name + '.mat']
    return mdict


def write_mat(out_to, data):
    sio.savemat(out_to, mdict=dict(data=data))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('path', type=str, help='path to a folder containing the .mat files to merge')
    p.add_argument('--output', '-o', type=str, default=None,
                   help='output file name, default is the name of the input folder')
    p.add_argument('var_regexp', type=str,
                   help=('Name matching python regexps that also defines the named variables. Semicolon-separated, '
                         'first matching regexp counts. '
                         'e.g. (?P<mod>\w+)_(?P<time>\d+)\.bw;(?P<mod>[\w-]+)_(?P<time>\d)\.bw'))
    args = p.parse_args()
    args.__dict__['var_regexp'] = [re.compile(x) for x in args.var_regexp.split(';')]
    if args.output is None:
        args.__dict__['output'] = os.path.split(args.path)[1] + '.mat'
    return args

if __name__ == '__main__':
    args = parse_args()
    vars, file_map, file_list, var_lists = get_data(args.path, args.var_regexp)
    mdict = build_matlab_objects(vars, file_map, file_list, var_lists)
    print mdict