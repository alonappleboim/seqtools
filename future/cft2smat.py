import sys
import scipy.io as sio
import numpy as np
from scipy import sparse
import argparse


INTERPRETER = '/cs/bd/tools/nflab_env/bin/python3.4'
if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    import subprocess as sp
    import os
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()


def parse_chrlen(chrlen_file):
    cl = {}
    with open(chrlen_file) as cl_file:
        for line in cl_file:
            c, l = line.strip().split('\t')
            cl[c] = int(l)
    return cl


def read_cft(cft_in):
    data = {}
    for line in cft_in:
        chr, fr, to = line.strip().split('\t')
        if chr not in data:
            data[chr] = ([],[])
        data[chr][0].append(int(fr))
        data[chr][1].append(int(to))
    return data


def write_mat(out_to, data, strand, meta, chr_lengths):
    s = {}
    for chr, (frs, tos) in data.items():
        s[chr] = sparse.csc_matrix((np.ones(len(frs)), (np.asarray(frs)-1, np.asarray(tos)-1)),
                                   shape=(chr_lengths[chr],chr_lengths[chr]), dtype=float)
    s['meta'] = meta
    s['strand'] = strand
    sio.savemat(out_to, mdict=dict(data=s))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input', '-i', type=str, default=None,
                   help='cft input file: chr\tfrom\tto, default is stdin')
    p.add_argument('name', type=str, help='sample name')
    p.add_argument('--output', '-o', type=str, default=None,
                   help='to which .mat output is written. default is name.cft.mat')
    p.add_argument('--strand', '-s', type=str, choices=['w','c','no'], default='no',
                   help='cft strand ("no" if not strand-specific)')
    p.add_argument('--chr_len_file', '-clf', type=str, default='/cs/wetlab/genomics/scer/genome/sacCer3.sizes',
                   help='path to chromosome lengths file')
    args = p.parse_args()
    if args.input is None:
        args.__dict__['input'] = sys.stdin
    else:
        args.__dict__['input'] = open(args.input)

    if args.output is None:
        args.__dict__['output'] = args.name + '.cft.mat'

    return args

if __name__ == '__main__':
    args = parse_args()
    write_mat(args.output, read_cft(args.input), args.strand, args.name, parse_chrlen(args.chr_len_file))
