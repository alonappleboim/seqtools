#! /cs/bd/tools/nflab_env/bin/python3.4

import sys
import os
import re
import shlex
import subprocess as sp
import scipy.io as sio
import numpy as np
from scipy import sparse
import argparse
from collections import OrderedDict

BW2W = '/cs/bd/tools/bigWigToWig'
if not os.path.exists(BW2W):
    BW2W = '/Users/user/nflab_scripts/Shell/bigWigToWig'

# sys.path.append(os.path.split(os.path.split(__file__)[0])[0])
# INTERPRETER = '/cs/bd/tools/nflab_env/bin/python3.4'
# if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
#     import subprocess as sp
#     import os
#     scriptpath = os.path.abspath(sys.modules[__name__].__file__)
#     sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
#     exit()


def cpv_iter(bw_in, chr_map, tmp_suff='.tmp'):
    """
    given a bigwig file, iterate over (chr, pos, val) triplets
    :param bw_in: bigwig input
    :param tmp_name: temporary file suffix
    :param chr_map: an OrderedDict of chr -> length
    :return: yields triplets of chromsome, position, value over all non-zero positions
    """
    def parse_fixed_step(input, chr_map):
        rm, csl, chr = 0, 0, None
        with open(input) as IN:
            for i, line in enumerate(IN):
                if line.startswith('fixedStep') or line.startswith('variableStep'):
                    rm += 1  # remove line from position count
                    c = line.strip().split(' ')[1].split('=')[1]
                    if c != chr:  # new chromsome
                        chr = c
                        l = chr_map[chr]
                        csl += l
                else:
                    val = float(line.strip())
                    if val > 0:
                        pos = i - rm - csl + l
                        yield chr, int(pos), float(val)

    def parse_variable_step(input, chr_map):
        rm, csl, chr = 0, 0, None
        with open(input) as IN:
            for line in IN:
                if line.startswith('variableStep'):
                    chr = line.strip().split(' ')[1].split('=')[1]
                else:
                    pos, val = line.strip().split('\t')
                    pos = int(pos)
                    yield chr, int(pos), float(val)

    def parse_bed_graph(input, chr_map):
        with open(input) as IN:
            for i, line in enumerate(IN):
                if line.startswith('#bedGraph'): continue
                chr, frm, to, val = line.strip().split('\t')
                for p in range(int(frm), int(to)):
                    yield chr, int(p), float(val)

    tmp_name = bw_in + tmp_suff
    comm = '%s %s %s' % (BW2W, bw_in, tmp_name)
    sp.Popen(shlex.split(comm)).communicate()
    t = open(tmp_name)
    hdr = t.readline()
    t.close()
    if hdr.startswith('fixedstep'): iter = parse_fixed_step
    elif hdr.startswith('variableStep'): iter = parse_variable_step
    elif hdr.startswith('#bedGraph'): iter = parse_bed_graph
    else: raise IOError('unknown Wig format: %s' % hdr)
    for cpv in iter(tmp_name, chr_map): yield cpv
    os.remove(tmp_name)


def organize_files(path, var_parsers):
    """
    Parse the files in path in light of given regexps

    :param path: a path containing files corresponding to given regexps
    :param var_parsers: a list of regexps that capture variables
    :return:
        vars - an ordered dictionary of sample variables, mapping to ordered lists of variable values
        file_list - a list of files that matched criteria, in order of iteration
        var_lists - a map of lists, one per variable, of variable values, corresponding to files in file_list
    """
    file_list, var_lists, vars = [], {}, OrderedDict()
    for vp in var_parsers:
        for v in vp.groupindex.keys():
            if v not in vars:
                vars[v] = []
                var_lists[v] = []
    for f in os.listdir(os.path.abspath(path)):
        for vp in var_parsers:
            m = vp.match(f)
            if m is not None: break
        if m is None: continue  # no match
        fpath = path + os.path.sep + f
        s_dict, s_vals = m.groupdict(), []
        for v in vars:
            val = s_dict[v] if v in s_dict else None
            var_lists[v].append(val)
            s_vals.append(val)
        for var, val in s_dict.items():
            if val not in vars[var]: vars[var].append(val)
        file_list.append(fpath)
    return vars, file_list, var_lists


def bw2py(bw_in, chr_map):
    """
    Read a bw file and convert it to a map: chr -> (positions, values)
    :param bw_in: bigwig (.bw) file path
    :param chr_map:  A chr -> length map
    :return: map chr -> (positions, values)
    """
    map = {c: ([],[]) for c in chr_map}
    for c, p, v in cpv_iter(bw_in, chr_map):
        map[c][0].append(p)
        map[c][1].append(v)
    return map


def build_matlab_objects(vars, file_list, var_lists, args):
    """
    prep the output MATLAB objects

    :param vars: an OrderedDict of variable name -> variable values
    :param file_list: an ordered list of files
    :param var_lists: a corresponding list of variable values per file
    :return:
        mdict - a MATLAB struct with following fields:
                * d - a cell array of #chromosomes, each holding a matrix of N_ixS, where N_i is the chromosome
                      length, and S is the number of samples.
                * l - a legend for the data, a struct, with following fields:
                      * samples - a list of sample names in the order they appear in "d{i}"
                      * sample_vars - a struct with every field corresponding to a variable of the samples, mapping to a
                                      list of values of that variable in the order of samples in "d"
                      * vars - a struct with variables values in order
                      * c - chromosome names, in order of "d"
    """
    S, C = len(file_list), len(args.chr_map)
    # init = sparse.csc_matrix if as_sparse else np.zeros
    d = [np.zeros((L,S)) for L in args.chr_map.values()] if not args.sparse else [([],[],[]) for _ in args.chr_map]
    out = {'d': d,
           'l': {'samples': np.asarray([os.path.split(f)[1][:-3] for f in file_list], dtype='object'),
                 'sample_vars': {v: np.asarray(var_lists[v],dtype='object') for v in vars.keys()},
                 'vars': {v : np.asarray(vals, dtype='object') for v, vals in vars.items()},
                 'c': np.asarray(args.chr_map.keys(),dtype='object')}
           }
    for si, f in enumerate(file_list):
        if args.verbose:
            sys.stderr.write('processing %s...\n' % f)
        fdata = bw2py(f, args.chr_map)
        for chi, (ch, _) in enumerate(args.chr_map.iteritems()):
            if ch not in fdata: continue
            pos, vals = fdata[ch]
            if not args.sparse:
                out['d'][chi][pos,si] = vals
            else:
                out['d'][chi][0].extend(pos) #rows
                out['d'][chi][1].extend([si]*len(pos)) #cols
                out['d'][chi][2].extend(vals) #values
    if args.sparse:
        for chi, (chr, chL) in enumerate(args.chr_map.iteritems()):
            R, C, V = out['d'][chi]
            out['d'][chi] = sparse.csc_matrix((V, (R, C)), shape=(chL,S), dtype=float)
    return {args.outname: out}


def reorder_samples(file_list, var_lists, vars, args):
    for v in args.order_by[::-1]:
        ord = np.argsort(var_lists[v],kind='mergesort')
        var_lists = {v: [vl[i] for i in ord] for v, vl in var_lists.items()}
        file_list = [file_list[i] for i in ord]
    vars = OrderedDict([(v,sorted(vars[v])) for v in args.order_by])
    return file_list, var_lists, vars


def parse_args():

    def parse_chrmap(chr_map):
        cm = OrderedDict()
        if os.path.exists(chr_map):
            with open(chr_map) as CM:
                for line in CM:
                    c, l = line.strip().split('\t')
                    cm[c] = int(l)
        else:
            for cl in chr_map.split(';'):
                c, l = cl.split(':')
                cm[c] = int(l)
        return cm

    def build_default_regexp(path):
        vars = set([])
        pat_str = ''
        for f in os.listdir(path):
            if not f.endswith('.bw'): continue
            for var_val in f.replace('.bw','').split('_'):
                if '-' not in var_val: continue
                var, _ = var_val.split('-')
                if var in vars: continue
                vars.add(var)
                pat_str += '%s-(?P<%s>[\w\.]+)_' % (var, var)
        return re.compile(pat_str[:-1] + '.bw'), list(vars)

    p = argparse.ArgumentParser()
    p.add_argument('path', type=str, help='path to a folder containing the bw files to merge')
    p.add_argument('outname', type=str, help=('The name of the resulting matlab struct'))
    p.add_argument('--verbose', '-v', action='store_false', help='do not output processing info to stderr')
    p.add_argument('--sparse', '-ns', action='store_false', help='should output NOT be in sparse format')
    p.add_argument('--output', '-o', type=str, default=None,
                   help='output file name, default is stdout')
    p.add_argument('--var_regexp', '-ve', type=str, default=None,
                   help=('Name matching python regexps that also defines the named variables. Semicolon-separated, '
                         'first matching regexp counts. '
                         'e.g. (?P<mod>\w+)_(?P<time>\d+)\.bw;(?P<mod>[\w-]+)_(?P<time>\d)\.bw '
                         'If not given, the assumed structure is <varname>-<varval>(_<varname>-<varval>)*bw'))
    p.add_argument('--order_by', '-ob', type=str, default=None,
                   help=('comm-separated list of variables (as in var_regexp)'))
    p.add_argument('--chr_map', '-cm', type=str, default='/cs/wetlab/genomics/scer/genome/sacCer3_ordered.sizes',
                   help=('A list of chr1_name:chr1_length;chr2_name:chr2_length... for the chromosomes in input bw '
                         'files, the given order will determine the order in the output cell array. If a valid file '
                         'path is given, the data is read from file, one line per chr, tab delimited name\tlength.'))
    args = p.parse_args()
    if args.output is None:
        args.__dict__['output'] = sys.stdout
    else:
        args.__dict__['output'] = open(args.__dict__['output'], 'wb')
    args.__dict__['chr_map'] = parse_chrmap(args.chr_map)
    if args.order_by: args.__dict__['order_by'] = args.order_by.split(',')
    if not args.var_regexp or args.var_regexp is None:
        ve, vord = build_default_regexp(args.path)
        args.__dict__['var_regexp'] = [ve]
        if args.order_by is None: args.__dict__['order_by'] = vord
    else:
        args.__dict__['var_regexp'] = [re.compile(x) for x in args.var_regexp.split(';')]
        if args.order_by is None:
            vars = set([])
            for ve in args.var_regexp: vars |= set(ve.groupindex.keys())
            args.__dict__['order_by'] = list(vars)
    return args


if __name__ == '__main__':
    args = parse_args()
    vars, file_list, var_lists = organize_files(args.path, args.var_regexp)
    file_list, var_lists, vars = reorder_samples(file_list, var_lists, vars, args)
    if args.verbose:
        sys.stderr.write('marging %i files...\n' % len(file_list))
    mdict = build_matlab_objects(vars, file_list, var_lists, args)
    sio.savemat(args.output, mdict)

