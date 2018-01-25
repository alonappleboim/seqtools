#! /cs/bd/tools/nflab_env/bin/python3.4

# collect all experiment statistics, load mat Ts4tU.mat and cov4tU.mat, and consolidate everything into a single MATLAB
# struct.

import sys
import os
import re
import pickle
import scipy.io as sio
import numpy as np
from scipy import sparse
import argparse
from collections import OrderedDict


def isint(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


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
    sys.stderr.write('Processing %i index files.\n' % len(file_list))
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

    def sum_dups(file):
        s = 0
        try:
            with open(file) as F:
                for l in F:
                    u, c = [int(x) for x in l.strip().split('\t')]
                    s += (u - 1) * c
        except IOError as e:
            sys.stderr.write(str(e) + '\n')
        return s

    def parse_align_stats(file, buffer):
        try:
            with open(file) as F:
                for l in F:
                    ''' example of a file
                    1194660 reads; of these:
                    1194660 (100.00%) were unpaired; of these:
                    246134 (20.60%) aligned 0 times
                    742514 (62.15%) aligned exactly 1 time
                    206012 (17.24%) aligned >1 times
                    79.40% overall alignment rate
                    '''
                    l = l.strip()
                    if "reads; of these" in l: buffer[0] = float(l.split(' ')[0])
                    if "aligned 0 times" in l: buffer[1] = float(l.split(' ')[0])
                    if "aligned exactly" in l: buffer[2] = float(l.split(' ')[0])
                    if "aligned >1" in l: buffer[3] = float(l.split(' ')[0])
                    if 'overall' in l: buffer[4] = float(l.split(' ')[0].replace('%', ''))
        except (IOError, ValueError) as e:
            sys.stderr.write(str(e) + '\n')

    nF = len(files['nrm-align'])
    p_align = np.zeros((nF, 5))  # preliminary alignment, just for total reads
    s_align = np.zeros((nF, 5))  # spike-in alignment: total, not aligned, unique, multiple, overall
    n_align = np.zeros((nF, 6))  # normal alignment: total, not aligned, unique, multiple, overall, dup
    t_align = np.zeros((nF, 6))  # 3-letter alignment: total, not aligned, unique, multiple, overall, dup
    rcu = np.zeros((nF, 3))
    rcu_ord = ['rejected', 'converted', 'unconverted']
    for i, f in enumerate(files['prelim-align']): parse_align_stats(f, p_align[i, :])
    for i, f in enumerate(files['spk-align']): parse_align_stats(f, s_align[i, :])
    for i, f in enumerate(files['nrm-align']): parse_align_stats(f, n_align[i,:])
    for i, f in enumerate(files['nrm-dhist']): n_align[i, 5] = sum_dups(f)
    for i, f in enumerate(files['tlg-align']): parse_align_stats(f, t_align[i, :])
    for i, f in enumerate(files['tlg-dhist']): t_align[i, 5] = sum_dups(f)
    for i, f in enumerate(files['cl']):
        try:
            with open(f) as F:
                for l in F:
                    n, v = l.split('\t')
                    j = rcu_ord.index(n)
                    rcu[i,j] = float(v)
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
    d = OrderedDict()
    # taking only "total" reads since this is only a preliminary
    # alignment to filter out spike-in aligned reads:
    d['total'] = p_align[:, 0]
    # taking "not aligned" reads since this is alignment of spike-in aligned
    # reads to s.cer, i.e. we want only reads that DO NOT align to cerevisiae:
    d['spk_only'] = s_align[:, 1]
    a_ord = ['not_aligned', 'unique', 'multiple', 'overall', 'dup']
    for i, name in enumerate(a_ord):
        d['align_' + name] = n_align[:, i + 1] # i+1 -> skip total, which is first in align stats
        d['cnv_align_' + name] = t_align[:, i + 1]
    for i, name in enumerate(rcu_ord): d[name] = rcu[:, i]
    return d


def merge_t_stats(files):
    d = []
    for i, f in enumerate(files):
        try:
            d.append(np.loadtxt(f, skiprows=1, dtype=float))
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
            d.append(np.zeros((0, 0)))
    lg = np.asarray(['observed', 'converted', 'count'], dtype='object')
    return lg, np.array(d, dtype='object')


def merge_dup_stats(files):
    d = []
    for i, f in enumerate(files):
        try:
            d.append(np.loadtxt(f, dtype=float))
        except IOError as e:
            sys.stderr.write(str(e) + '\n')
            d.append(np.array([]))
    lg = np.asarray(['counts', 'unique molecules'], dtype='object')
    return lg, d


def merge_annotation_counts(files, annot_file):
    lg = {}
    try:
        # parse tts file for legend
        afile = open(annot_file)
        hdr = afile.readline().strip().split(' ')
        lg['W'] = np.array([float(hdr[1]), float(hdr[2])])
        lg['acc'] = np.asarray([line.strip().split('\t')[3] for line in afile], dtype='object')
        N = len(lg['acc'])

        # and iterate over all count files to merge them
        d = np.zeros((N, len(files), 3))
        for fi, f in enumerate(files):
            fd = np.loadtxt(f, dtype=int, skiprows=1)
            for ti in range(3):
                d[fd[:,0]-1, fi, ti] = fd[:,ti+1]
        lg['data_types'] = np.asarray(open(f).readline().strip().split('\t')[1:], dtype='object')
    except IOError as e:
        sys.stderr.write(str(e) + '\n')
        d = np.zeros((0, 0))
    return {'l':lg, 'd':d}


def merge_mut_stats(files):
    lg = {'alphabet': np.asarray(['A','C','G','T','N'],dtype='object'),
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


def merge_ahists(files):
    lg = np.asarray(['#A','cnt'],dtype='object')
    ds, mxd = [], 0
    for i, f in enumerate(files):
        try:
            with open(f) as F:
                fd = np.loadtxt(f, dtype=int)
                ds.append(fd)
                mxd = max(mxd, fd[:,0].max())
        except IOError as e:
            sys.stderr.write(str(e)+'\n')
    d = np.zeros((mxd, len(files)))
    for i, di in enumerate(ds):
        d[di[:,0]-1, i] = di[:,1]
    return lg, d


def collect_stats(files, var_lists, vars):
    stats = {}
    stats['general'] = merge_general_stats(files)
    stats['read_t'] = {}
    stats['read_t']['l'], stats['read_t']['d'] = merge_t_stats(files['rt'])
    stats['genomic_t'] = {}
    stats['genomic_t']['l'], stats['genomic_t']['d'] = merge_t_stats(files['gt'])
    stats['polya_hist'] = {}
    stats['polya_hist']['l'], stats['polya_hist']['d'] = merge_ahists(files['ah'])
    stats['ndup'] = {}
    stats['ndup']['l'], stats['ndup']['d'] = merge_dup_stats(files['nrm-dhist'])
    stats['dup'] = {}
    stats['dup']['l'], stats['dup']['d'] = merge_dup_stats(files['tlg-dhist'])
    stats['snp'] = {}
    stats['snp']['l'], stats['snp']['d'] = merge_mut_stats(files['mut'])

    stats['l'] = {'samples': np.asarray([os.path.split(f)[1][:-len('.index')] for f in files['index']], dtype='object'),
                  'sample_vars': {v: np.asarray(var_lists[v],dtype='object') for v in vars.keys()},
                  'vars': {v : np.asarray(vals, dtype='object') for v, vals in vars.items()}}
    return stats


def collect_pconv_fit(files):
    d = np.stack([np.loadtxt(f) for f in files], axis=2)
    fit = pickle.loads(open(f).readline().replace('#','').strip().replace(';','\n'))
    fit['d'] = d
    return fit


def fit_params(files, stats):
    fit = {}
    fit['p_conv'] = collect_pconv_fit(files)
    return fit


def file_map(file_list):
    fmap = {'index':[]}
    for IF in file_list:
        fmap['index'].append(IF)
        for line in open(IF):
            ftype, fpath = line.rstrip().split(':')
            if ftype not in fmap: fmap[ftype] = []
            fmap[ftype].append(fpath)
    return fmap


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
    p.add_argument('--output', '-o', type=str, default=None, help='output file name, default is stdout')
    p.add_argument('--annot_file', '-af', type=str, default='STATS/tts_w.bed',
                   help='the bed file used to count reads into annotations.')
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
    fmap = file_map(file_list)
    stats = collect_stats(fmap, var_lists, vars)
    fit = fit_params(fmap['tbino_fit'], stats)
    acnt = {} if 'annot_cnt' not in fmap else merge_annotation_counts(fmap['annot_cnt'], args.annot_file)
    sio.savemat(args.output, {args.outname: dict(fit=fit, stats=stats, acnt=acnt)})

