"""
Reads a [chr,pos,strand] file (from stdin) and outputs read statistics per entry according to
specified format in command line arguments.
"""

import argparse
import os
import sys
from abc import ABCMeta, abstractmethod

project_dir = os.path.sep.join(sys.modules[__name__].__file__.split(os.path.sep)[:-2])
sys.path.append(project_dir)

from common.config import *

if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    import subprocess as sp
    import os
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()

import numpy as np
import pysam
import scipy.io as sio

from common.utils import chr_lengths

MAX_INSERT = 1001


class Segment(object):

    def __init__(self, r, should_pair=False, afile=None):
        if should_pair:
            rr = afile.mate(r)
            self.fr = min(r.reference_start, rr.reference_start)
            self.to = max(r.reference_end, rr.reference_end)
        else:
            self.fr = r.reference_start
            self.to = r.reference_end

    def __len__(self):
        return self.to - self.fr


class SegmentTransform(object):
    __metaclass__ = ABCMeta

    name = None
    args = {}
    desc = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @abstractmethod
    def transform(self, s): pass

    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)


class Step(SegmentTransform):
    name = 'step'
    args = {'w': (int, 1, 'half-width to extend from center (to each side)')}
    desc = 'step function around read center'

    def transform(self, s):
        l = len(s)
        t = np.zeros(l)
        c = round(l/2)
        t[max(0,c-self.w):min(c+self.w,l)] = 1
        return t


class Coverage(SegmentTransform):
    name = 'cov'
    args = {}
    desc = 'uniform 1 coverage'

    def transform(self, s):
        return np.ones(len(s))


class FivePrime(SegmentTransform):
    name = '5p'
    args = {}
    desc = "indicator at read start (5')"

    def transform(self, s):
        t = np.zeros(len(s))
        t[0] = 1
        return t


class ThreePrime(SegmentTransform):
    name = '3p'
    args = {}
    desc = "indicator at read end (3')"

    def transform(self, s):
        t = np.zeros(len(s))
        t[-1] = 1
        return t


class Spike(SegmentTransform):
    name = 'spike'
    args = {'a': (float, 5, 'the steepness of the pike, the higher the more spiky the transform')}
    desc = "a concave spike at the read center"

    def __init__(self, **kwargs):
        super(Spike, self).__init__(**kwargs)
        self.halfspike = 1 / (np.array(range(1,round(MAX_INSERT/2)+1)) ** self.a)

    def transform(self, s):
        sl = len(s)
        h1 = self.halfspike[:round(sl/2)]
        if sl % 2 == 0:
            t = np.hstack((h1[::-1], h1))
        else:
            t = np.hstack((h1[::-1], h1[1], h1))
        return t


def reshape(bam_path, annot_file, transform, out_name, output_file, win=[-500,250], is_paired=True, same_strand=True):
    bam_in = pysam.AlignmentFile(bam_path)
    chrlens = chr_lengths()
    vectors = []
    annot = []
    for line in annot_file:
        id, chr, pos, strand = line.strip().split(ANNOT_DELIM)
        annot.append(id)
        pos, strand = int(pos), (-1) ** (1 - (strand == '+'))
        is_rev = strand == -1
        fr, to = [pos + strand*x for x in win[::strand]]
        fd = abs(min(0, fr))  # >0 only in case that reached end of chromosome
        fr, to = max(0, fr), min(to, chrlens[chr])
        w = np.zeros(win[1] - win[0]+1)
        wl = len(w)
        for r in bam_in.fetch(chr, fr, to):
            if is_paired:
                if not r.is_read1 or not r.is_proper_pair: continue
                s = Segment(r, should_pair=True, afile=bam_in)
            else:
                if not same_strand and r.is_reverse == is_rev: continue
                if same_strand and not r.is_reverse == is_rev: continue
                s = Segment(r)
            add = transform(s)
            if is_rev: add = add = add[::-1]
            if strand == 1:
                wfr, wto = s.fr - fr, s.to - fr
            else:
                wfr, wto = to - s.to, to - s.fr
            sfr, sto = abs(min(wfr, 0)), len(s) - abs(min(len(w) - wto, 0))
            wfr, wto = max(0, wfr), min(wto, wl)
            w[fd+wfr:fd+wto] += add[sfr:sto]
        vectors.append(w)
    s = {'d': np.vstack(vectors), 'w': w, 'l': annot}
    sio.savemat(output_file, {out_name: s})


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
    g.add_argument('--is_paired', '-p', action='store_true',
                   help='whether data is paired reads or not.')
    g.add_argument('--not_same_strand', '-s', action='store_true',
                   help='should directionality be flipped relative to annotation?')
    g.add_argument('--annot_in', '-ain', type=str, default=None,
                   help='path to a "id\tchr\tpos\tstrand" annotation file, default is stdin')
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
            if issubclass(transform_module.__dict__[x], SegmentTransform):
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
    args = parse_arguments(build_parser())
    print('reshaping %s...' % args.bam_in)
    reshape(args.bam_in, args.annot_in, args.transform,
            args.out_name, args.output_file, args.w, args.is_paired, not args.not_same_strand)
    print('wrote file: %s' % args.output_file)