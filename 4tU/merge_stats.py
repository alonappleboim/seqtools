#! /cs/bd/tools/nflab_env/bin/python3.4

# collect all experiment statistics, load mat Ts4tU.mat and cov4tU.mat, and consolidate everything into a single MATLAB
# struct.

import sys
import os
import re
import scipy.io as sio
import numpy as np
from scipy import sparse
import argparse
from collections import OrderedDict


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
    for v, type in args.order_by[::-1]:
        vals = [float(vv) for vv in var_lists[v]] if type == 'num' else var_lists[v]
        ord = np.argsort(vals,kind='mergesort')
        var_lists = {v: [vl[i] for i in ord] for v, vl in var_lists.items()}
        file_list = [file_list[i] for i in ord]
    vars = OrderedDict([(v,sorted(vars[v])) for v, _ in args.order_by])
    return file_list, var_lists, vars


def merge_general_stats(files):
    lg = ['total', 'nrm-aligned', 'rejected', 'converted', 'unconverted', 'nrm-dup', 'dup']
    total, naln, rej, cnv, unc, ndup, dup = [], [], [], [], [], [], []
    for f in files['nrm-align']:
        try:
            with open(f) as F:
                for i, l in enumerate(F):
                    if i == 0: total.append(int(l.strip().split(' ')[0]))
                    if 'overall' in l: naln.append(float(l.split(' ')[0].replace('%','')))
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
            total.append(0)
            naln.append(0)
    for f in files['cl']:
        try:
            with open(f) as F:
                for l in F:
                    if l.startswith('converted'): cnv.append(int(l.split('\t')[1]))
                    if l.startswith('unconverted'): unc.append(int(l.split('\t')[1]))
                    if l.startswith('rejected'): rej.append(int(l.split('\t')[1]))
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
            cnv.append(0)
            unc.append(0)
            rej.append(0)
    for f in files['tlg-dhist']:
        c = 0
        try:
            with open(f) as F:
                for l in F:
                    if l.startswith('1'): continue
                    c += int(l.split('\t')[1])
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
        dup.append(c)
    for f in files['nrm-dhist']:
        c = 0
        try:
            with open(f) as F:
                for l in F:
                    if l.startswith('1'): continue
                    c += int(l.split('\t')[1])
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
        ndup.append(c)
    d = np.asarray([total, naln, rej, cnv, unc, ndup, dup]).T
    return np.asarray(lg,dtype='object'), d


def merge_t_stats(files):
    lg = ['dim1=observed', 'dim2=converted', 'dim3=ordered-samples']
    samples = []
    mxo, mxc = 0, 0
    for i, f in enumerate(files):
        try:
            d = np.loadtxt(f, skiprows=1, dtype=int)
            if d.shape[0] == 0: raise IOError("%s is an empty file" % f)
            if len(d.shape) == 1: d = d.reshape([1, 3])
            mxo = max(mxo, max(d[:, 0]))
            mxc = max(mxc, max(d[:, 1]))
            samples.append(d)
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
            samples.append(None)
    d = np.zeros([mxo+1,mxc+1,i+1], dtype=float)
    for i, s in enumerate(samples):
        if s is None: continue
        si = i*np.ones(s.shape[0],dtype=int)
        d[s[:, 0], s[:, 1], si] = s[:, 2].astype(float)
    return np.asarray(lg,dtype='object'), d


def merge_dup_stats(files):
    lg = ['dim2=dup-hist','dim1=ordered-samples']
    samples = []
    mxd = 0
    for i, f in enumerate(files):
        try:
            d = np.loadtxt(f, dtype=int)
            if len(d.shape) == 1: d = d.reshape([1, 2])
            samples.append(d)
            mxd = max(mxd, max(d[:, 0]))
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
            samples.append(None)
    d = np.zeros([mxd+1,i+1], dtype=float)
    for i, s in enumerate(samples):
        if s is None: continue
        d[s[:, 0], i] = s[:, 1].astype(float)
    return np.asarray(lg,dtype='object'), d


def merge_mut_stats(files):
    lg = {'alphabt': np.asarray(['A','C','G','T','N'],dtype='object'),
          'dims': np.asarray(['ref','obs','strand','samples'], dtype='object'),
          'strand': np.asarray(['W','C'], dtype='object')}
    amap = {'A':0,'C':1,'G':2,'T':3,'N':4}
    smap = {'W':0,'C':1}
    d = np.zeros([len(amap),len(amap),len(smap),len(files)])
    for i, f in enumerate(files):
        try:
            with open(f) as F:
                F.readline() # get rid of header
                for l in F:
                    s, r, o, c = l.strip().split('\t')
                    d[amap[r],amap[o],smap[s],i] = float(c)
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
    return lg, d


def collect_stats(index_files, var_lists, vars):
    files = {}
    for IF in index_files:
        for line in open(IF):
            ftype, fpath = line.rstrip().split(':')
            if ftype not in files: files[ftype] = []
            files[ftype].append(fpath)
    stats = {}
    stats['general'] = {}
    stats['general']['l'], stats['general']['d'] = merge_general_stats(files)
    stats['read_t'] = {}
    stats['read_t']['l'], stats['read_t']['d'] = merge_t_stats(files['rt'])
    stats['genomic_t'] = {}
    stats['genomic_t']['l'], stats['genomic_t']['d'] = merge_t_stats(files['gt'])
    stats['ndup'] = {}
    stats['ndup']['l'], stats['ndup']['d'] = merge_dup_stats(files['nrm-dhist'])
    stats['dup'] = {}
    stats['dup']['l'], stats['dup']['d'] = merge_dup_stats(files['tlg-dhist'])
    stats['snp'] = {}
    stats['snp']['l'], stats['snp']['d'] = merge_mut_stats(files['mut'])

    stats['l'] = {'samples': np.asarray([os.path.split(f)[1][:-len('.index')] for f in index_files], dtype='object'),
                  'sample_vars': {v: np.asarray(var_lists[v],dtype='object') for v in vars.keys()},
                  'vars': {v : np.asarray(vals, dtype='object') for v, vals in vars.items()}}

    return stats


def parse_args():

    def build_default_regexp(path):
        vars = set([])
        pat_str = ''
        for f in os.listdir(path):
            if not f.endswith('.index'): continue
            for var_val in f.replace('.index','').split('_'):
                if '-' not in var_val: continue
                var, _ = var_val.split('-')
                if var in vars: continue
                vars.add(var)
                pat_str += '%s-(?P<%s>[\w\.]+)_' % (var, var)
        return re.compile(pat_str[:-1] + '.index'), list(vars)

    p = argparse.ArgumentParser()
    p.add_argument('--stats_dir', '-sd', type=str, help='path to statistics directory', default='STATS')
    p.add_argument('--cov_mat', '-cm', type=str, help='path to Matlab coverage matrix', default='tmp/cov.mat')
    p.add_argument('--t_mat', '-tm', type=str, help='path to Matlab genomic T statistics matrix', default='tmp/t.mat')
    p.add_argument('--output', '-o', type=str, default=None, help='output file name, default is stdout')
    p.add_argument('--outname', '-on', type=str, default='tU', help='output struct name, default is tU')
    p.add_argument('--var_regexp', '-ve', type=str, default=None,
                   help=('Name matching python regexps that also defines the named variables. Semicolon-separated, '
                         'first matching regexp counts. '
                         'e.g. (?P<mod>\w+)_(?P<time>\d+)\.bw;(?P<mod>[\w-]+)_(?P<time>\d)\.bw '
                         'If not given, the assumed structure is <varname>-<varval>(_<varname>-<varval>)*\.index'))
    p.add_argument('--order_by', '-ob', type=str, default=None,
                   help=('A list of variarbles and how to sort them. Order of variables determines order in resulting'
                         ' hub. e.g. "mod:str;time:num" will be first sorted by "mod"(string), then by "time"(number)'))
    args = p.parse_args()
    if args.output is None:
        args.__dict__['output'] = sys.stdout
    else:
        args.__dict__['output'] = open(args.output, 'wb')

    if args.order_by: args.__dict__['order_by'] = [x.split(':') for x in args.order_by.split(';')]

    if not args.var_regexp or args.var_regexp is None:
        ve, vord = build_default_regexp(args.stats_dir)
        args.__dict__['var_regexp'] = [ve]
        if args.order_by is None: args.__dict__['order_by'] = [(v, 'str') for v in vord]
    else:
        args.__dict__['var_regexp'] = [re.compile(x) for x in args.var_regexp.split(';')]
        if args.order_by is None:
            vars = set([])
            for ve in args.var_regexp: vars |= set(ve.groupindex.keys())
            args.__dict__['order_by'] = [(v, 'str') for v in vars]
    return args

    return args


if __name__ == '__main__':
    args = parse_args()
    vars, file_list, var_lists = organize_files(args.stats_dir, args.var_regexp)
    file_list, var_lists, vars = reorder_samples(file_list, var_lists, vars, args)
    stats = collect_stats(file_list, var_lists, vars)
    t_mat = sio.loadmat(args.t_mat)
    cov_mat = sio.loadmat(args.cov_mat)
    mdict = {args.outname:{'cov' : cov_mat['cov4tU'], 't' : t_mat['Ts4tU'], 'stats' : stats}}
    sio.savemat(args.output, mdict)

