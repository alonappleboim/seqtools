#! /cs/bd/tools/nflab_env/bin/python3.4

import sys
import os
import heapq as hq
import scipy.io as sio
import numpy as np
import argparse
import re
from collections import OrderedDict

# sys.path.append(os.path.split(os.path.split(__file__)[0])[0])
# INTERPRETER = '/cs/bd/tools/nflab_env/bin/python3.4'
# if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
#     import subprocess as sp
#     import os
#     scriptpath = os.path.abspath(sys.modules[__name__].__file__)
#     sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
#     exit()

MIN_DENSITY = 1.01


def build_data_map(matpath):
    """
    :param matpath:
    :return:
    D, an OrderedDict: sample -> chr -> data vector
    mR, a map sample -> mean reads
    """
    D = OrderedDict()
    data = sio.loadmat(matpath)
    for k in data:
        if not k.startswith('__'): break
    data = data[k]
    leg, data = data['l'], data['d']
    data = data[0][0][0]
    samples = leg[0][0][0]['samples'][0][0]
    chrs = leg[0][0][0]['c'][0][0]
    L, mR = 0, {s[0]: 0 for s in samples}
    for c, d in zip(chrs, data):
        L += d.shape[0]
        for i, s in enumerate(samples):
            if s[0] not in D: D[s[0]] = OrderedDict()
            D[s[0]][c[0]] = d[:,i]
            mR[s[0]] += abs(np.sum(d[:,i]))
    mR = {s: n/L for s, n in mR.iteritems()}
    return D, mR


def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indices of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx


def call_pAs(d, min_density, min_width, max_width):
    h = []
    idx = contiguous_regions(np.asarray(d.todense()).flatten()!=0)
    for i, (start, end) in enumerate(zip(idx[:,0], idx[:,1])):
        if end - start >= max_width: continue
        k = np.sum(d[start:end])/(end-start)
        hq.heappush(h, (1.0/k, i, i)) # min heap, so inverse density is used to pop the maximal density

    used = set([])
    while h:
        k, fi, ti = hq.heappop(h)
        pos = fr, to = idx[fi,0], idx[ti,1]
        if pos in used:
            used.remove(pos)
            continue

        w = to-fr
        ld = np.sum(d[idx[fi-1,0]:to]) / (to-idx[fi-1,0]) if fi > 0 else -1
        rd = np.sum(d[fr:idx[ti+1,1]]) / (idx[ti+1,1]-fr) if ti < idx.shape[0]-1 else -1
        pobj = None
        if ld > rd:
            if ld > min_density:
                used.add((idx[fi-1,0], idx[fi-1,1]))
                pobj = 1.0/ld, fi-1, ti
        else:
            if rd > min_density:
                used.add((idx[ti+1,0], idx[ti+1,1]))
                pobj = 1.0/rd, fi, ti+1
        if pobj is None: # extension didn't work
            if (1.0/k > min_density) and (w >= min_width) and (w <= max_width):
                yield fr, to, 1/k
        else:
            hq.heappush(h, pobj)


def call_and_output(D, mR, args):
    """

    :param D: an OrderedDict: sample -> chr -> data vector
    :param mR: sample -> mean reads
    :param args: perocess input arguemnts
    :return:
    if a matlab output is requested a dictionary with:
        d - density data
        c - chromosome
        s - start
        e - end

    """
    table = {}
    for si, (s, sdict) in enumerate(D.items()):
        if args.debug:
            if si >= args.debug['nsample']: continue
        args.output.write('@%s\n' % s)
        s_md = max(mR[s] * args.min_density_factor, MIN_DENSITY)
        for ci, (c, d) in enumerate(sdict.items()):
            if args.debug:
                if ci >= args.debug['nchr']: continue
            if args.debug:
                d = d[:min(args.debug['nbp'], d.shape[0])]
            if args.verbose:
                sys.stderr.write('processing %s, %s (density threshold: %.2f)\n' % (s, c ,s_md))
            for fr, to, v in call_pAs(np.abs(d), s_md, args.min_width, args.max_width):
                args.output.write('%s\t%i\t%i\t%.2f\n' % (c, fr, to + 1, v))
                if args.mat is not None:
                    k = ci+1, fr, to
                    if k not in table: table[k] = [0] * len(D)
                    table[k][si] = v
    pA = {}
    if args.mat is not None:
        skeys = sorted(table.keys(), key=lambda x: x[2])
        skeys = sorted(skeys, key=lambda x: x[1])
        skeys = sorted(skeys, key=lambda x: x[0])
        N = len(skeys)
        pA = {'c': np.array([ci for ci, _, _ in skeys], dtype=np.double).reshape((N,1)),
              's': np.array([s for _, s, _ in skeys], dtype=np.double).reshape((N,1)),
              'e': np.array([e for _, _, e in skeys], dtype=np.double).reshape((N,1)),
              'd': np.array([table[k] for k in skeys], dtype=np.double)}

    return pA


def parse_args():

    def parse_chrmap(chr_map):
        cm = OrderedDict()
        if os.path.exists(args.chr_map):
            with open(args.chr_map) as CM:
                for line in CM:
                    c, l = line.strip().split('\t')
                    cm[c] = int(l)
        else:
            for cl in args.chr_map.split(';'):
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
                pat_str += '%s-(?P<%s>\w+)_' % (var, var)
        return re.compile(pat_str[:-1] + '.bw')

    p = argparse.ArgumentParser()
    p.add_argument('input', type=str, help='a .mat file as the one generated by bwhub2mat')
    p.add_argument('--output', '-o', default=None,
                   help=('''to which output, in concatnated bed format is writte. default is stdout each sample is '
                         'separated in output by a @<sample-name\\breakline> header. This can be separated to'
                         'individual files with the following awk commnand: 
                         awk '{if($0~"@") {bfile=substr($0,2,length($0))".bed";} else {print > bfile;}}' '''))
    p.add_argument('--mat', '-m', type=str, default=None,
                   help=('a comma spearated pair: <path>,<name>. If given, a matlab struct with name is saved in '
                         'given path'))
    p.add_argument('--verbose', '-v', action='store_false', help='do not output processing info to stderr')
    p.add_argument('--debug', '-d', default=None,
                   help=('debug mode, a comma separated triplet: <#samples, #chromosomes, #bp>'))
    p.add_argument('--min_density_factor', '-mdf', type=float, default=3500,
                   help=('minimum density factor over the mean value across the genome, for the segment to be '
                         'considered a pA site'))
    p.add_argument('--min_width', '-nw', type=int, default=1,
                   help='minimum width to be considered as a pA site')
    p.add_argument('--max_width', '-xw', type=int, default=100,
                   help='minimum width to be considered as a pA site')
    args = p.parse_args()

    if args.debug is not None:
        ns, nc, nb = [int(x) for x in args.debug.split(',')]
        args.__dict__['debug'] = {'nchr': nc, 'nsample': ns, 'nbp': nb}

    args.__dict__['output'] = sys.stdout if args.output is None else open(args.output, 'w')
    if args.mat is not None:
        args.__dict__['mat'] = args.mat.split(',')
    return args


def consolidate_segments(pA):
    """
    Align a collection of given segments to prevent minor differences in border definitions between samples
    e.g. chr IV: 100300-100310 + 100302-100311 => chr IV: 100300-100311.
    :param pA: as produced by the call_and_output function (sorted by chromosome, position)
    :return: pA, an updated consolidated map, with less segments, and updated densities
    """
    pA_u = {'c': [], 's':[], 'e':[], 'd': []}

    def add_segment(idx):
        s, e = np.min(pA['s'][idx]), np.max(pA['e'][idx])
        c = pA['c'][idx[0]]
        d = np.sum(pA['d'][idx,:] * np.tile(pA['e'][idx]-pA['s'][idx],(1,pA['d'].shape[1])),0) / (e - s)
        pA_u['c'].append(c)
        pA_u['s'].append(s)
        pA_u['e'].append(e)
        pA_u['d'].append(d)

    idx, i, s, e, c = [0], 1, pA['s'][0], pA['e'][0], pA['c'][0]
    while i < pA['c'].shape[0]:
        if pA['c'][i] != c or pA['s'][i] > e or pA['e'][i] < s:  # no overlap, segment finished
            add_segment(idx)
            idx, i, s, e, c = [i], i+1, pA['s'][i], pA['e'][i], pA['c'][i]  # init segment variables
        else:  # update segment variables
            s, e = min(s, pA['s'][i]), max(e, pA['e'][i])
            idx.append(i)
            i += 1

    N = len(pA_u['c'])
    return {'c': np.array(pA_u['c'], dtype=np.double).reshape((N, 1)),
            's': np.array(pA_u['s'], dtype=np.double).reshape((N, 1)),
            'e': np.array(pA_u['e'], dtype=np.double).reshape((N, 1)),
            'd': np.array(pA_u['d'], dtype=np.double)}

if __name__ == '__main__':
    args = parse_args()
    D, mR = build_data_map(args.input)
    pA = call_and_output(D, mR, args)
    pA = consolidate_segments(pA)
    if args.mat:
        sio.savemat(args.mat[0], mdict={args.mat[1]: pA})

