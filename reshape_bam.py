"""
Reads a [chr,pos,strand] file (from stdin) and outputs read statistics per entry according to
specified format in command line arguments.
"""

from abc import ABCMeta, abstractmethod
import pysam
import os
import argparse
import numpy as np
import scipy.io as sio
from config import *
from utils import chr_lengths

MAX_INSERT = 1001

class ReadTransform(object):
    __metaclass__ = ABCMeta

    name = None
    args = {}
    desc = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @abstractmethod
    def transform(self, r): pass

    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)


class Step(ReadTransform):
    name = 'step'
    args = {'w': (int, 2, 'half-width to extend from center (to each side)')}
    desc = 'step function around read center'

    def transform(self, r):
        t = np.zeros(r.reference_length)
        c = round(r.reference_length / 2)
        t[max(0,c-self.w):min(c+self.w,r.reference_length)] = 1
        return t


class Coverage(ReadTransform):
    name = 'cov'
    args = {}
    desc = 'uniform 1 coverage'

    def transform(self, r):
        return np.ones(r.reference_length)


class FivePrime(ReadTransform):
    name = '5p'
    args = {}
    desc = "indicator at read start (5')"

    def transform(self, r):
        t = np.zeros(r.reference_length)
        t[0] = 1
        return t


class ThreePrime(ReadTransform):
    name = '3p'
    args = {}
    desc = "indicator at read end (3')"

    def transform(self, r):
        t = np.zeros(r.reference_length)
        t[-1] = 1
        return t


class Spike(ReadTransform):
    name = 'spike'
    args = {'a': (float, 5, 'the steepness of the pike, the higher the more spiky the transform')}
    desc = "a concave spike at the read center"

    def __init__(self, **kwargs):
        super(Spike, self).__init__(**kwargs)
        self.halfspike = 1 / (np.array(range(1,round(MAX_INSERT/2)+1)) ** self.a)

    def transform(self, r):
        rl = r.reference_length
        h1 = self.halfspike[:round(rl/2)]
        if rl % 2 == 0:
            t = np.hstack((h1[::-1], h1))
        else:
            t = np.hstack((h1[::-1], h1[1], h1))
        return t


class Collector(object):

    def __init__(self, args):
        self.__dict__.update(args.__dict__)
        self.bam_in = pysam.AlignmentFile(args.bam_in)

    def collection(self):
        chrlens = chr_lengths()
        for line in self.annot_in:
            chr, pos, strand = line.strip().split(self.d)
            pos, strand = int(pos), (-1) ** (1 - (strand == '+'))
            is_rev = strand == -1
            fr, to = [pos + strand*x for x in self.w[::strand]]
            fd = abs(min(0, fr))  # >0 only in case that reached end of chromosome
            fr, to = max(0, fr), min(to, chrlens[chr])
            w = np.zeros(self.w[1] - self.w[0]+1)
            wl = len(w)
            for r in self.bam_in.fetch(chr, fr, to):
                if r.is_reverse == is_rev:
                    add = self.transform(r)
                    if is_rev: add = add = add[::-1]
                    if strand == 1:
                        wfr, wto = r.reference_start - fr, r.reference_end - fr
                    else:
                        wfr, wto = to - r.reference_end, to - r.reference_start
                    rfr, rto = abs(min(wfr, 0)), r.reference_length - abs(min(len(w) - wto, 0))
                    wfr, wto = max(0, wfr), min(wto, wl)
                    w[fd+wfr:fd+wto] += add[rfr:rto]
            yield w

    def collect(self):
        s = {'d': np.vstack(x for x in self.collection()),
             'w': self.w}
        sio.savemat(self.output_file, {self.out_name:s})


def parse_arguments(p):
    args = p.parse_args()
    if args.annot_in is None:
        args.__dict__['annot_in'] = sys.stdin
    else:
        args.__dict__['annot_in'] = open(args.annot_in)
    args.__dict__['w'] = [int(x) for x in args.w[1:-1].split(',')]

    if args.out_name is None:
        args.__dict__['out_name'] = os.path.split(args.output_file)[1].split(os.path.extsep)[0]

    args.__dict__['transform'] = transform_from_string(args.transform)
    return args


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('bam_in', type=str, default=None, help='path to an indexed bam file')
    g.add_argument('--annot_in', '-ain', type=str, default=None,
                   help='path to a "chr,pos,strand" annotation file, default is stdin')
    g.add_argument('--delim', '-d', type=str, default=',', dest='d',
                   help='annotation input delimiter')

    g = p.add_argument_group('Output')
    g.add_argument('--output_file', '-o', type=str, default='reshaped.mat',
                   help='to which output is written.')
    g.add_argument('--out_name', '-on', type=str, default=None,
                   help='name of matrix in matlab struct, default is the suffix of output_name (excluding .mat)')
    g.add_argument('--window', '-w', type=str, default='[-500,250]', dest='w',
                   help='Comma separated limits for the areas around annotations to be collected')
    g.add_argument('--transform', '-T', type=str, default='cov()',
                   help='the transform applied to each read. See -th for details.')

    return p


def collect_transforms():
    transforms = {}
    transform_module = sys.modules[__name__]
    for x in dir(transform_module):
        try:
            if issubclass(transform_module.__dict__[x], ReadTransform):
                tname = transform_module.__dict__[x].name
                if tname is not None:
                    transforms[tname] = transform_module.__dict__[x]
        except TypeError: pass
    return transforms


def transform_from_string(tstring):
    tclasses = collect_transforms()
    ts = []
    t, argdict = parse_transform_string(tstring)
    if t not in tclasses:
        raise ValueError('No transform %s' % t)
    T = tclasses[t]
    args = {}
    for aname, aval in argdict.items():
        args[aname] = T.args[aname][0](aval)
    for a, props in T.args.items():
        if a not in args: args[a] = props[1]  # default
    return T(**args)


def parse_transform_string(tstr):
    name, arg_s = tstr.split('(')
    arg_s = arg_s.strip()[:-1].split(',')
    args = {}
    for a in arg_s:
        a = a.strip()
        if not a: continue
        aname, aval = a.split('=')
        args[aname] = aval
    return (name, args)


if __name__ == '__main__':
    c = Collector(parse_arguments(build_parser()))
    c.collect()