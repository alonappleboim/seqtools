import sys
import argparse
import os
from common.config import *
import subprocess as sp
import shlex


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


def parse_fixedstep(input, out_f, cl):
    rm, csl, chr = 0, 0, None
    with open(input) as IN:
        for i, line in enumerate(IN):
            if line.startswith('fixedStep'):
                rm += 1 # remove line from position count
                c = line.strip().split(' ')[1].split('=')[1]
                if c != chr: # new chromsome
                    chr = c
                    l = cl[chr]
                    csl += l
            else:
                val = float(line.strip())
                if val > 0:
                    pos = i-rm-csl+l
                    out_f.write('%s\t%i\t%.2f\n' % (chr, pos, val))


def parse_bedgraph(input, out_f, cl):
    with open(input) as IN:
        for i, line in enumerate(IN):
            if line.startswith('#bedGraph'): continue
            chr, frm, to, val = line.strip().split('\t')
            for p in range(int(frm),int(to)):
                out_f.write('%s\t%i\t%s\n' % (chr, p, val))


def bw2bed(bw_in, cl, tmp_name, out_f):
    comm = '%s %s %s' % (EXEC['BW2W'], bw_in, tmp_name)
    sp.Popen(shlex.split(comm)).communicate()
    t = open(tmp_name)
    hdr = t.readline()
    t.close()
    if hdr.startswith('fixedstep'): parse_fixedstep(tmp_name, out_f, cl)
    elif hdr.startswith('#bedGraph'): parse_bedgraph(tmp_name, out_f, cl)
    os.remove(tmp_name)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('bw_in', type=str, help='BigWig input file')
    p.add_argument('--output', '-o', type=str, default=None,
                   help='to which bed output is written. default is stdout')
    p.add_argument('--chr_len_file', '-clf', type=str, default='/cs/wetlab/genomics/scer/genome/sacCer3.sizes',
                   help='path to chromosome lengths file')
    args = p.parse_args()
    if args.output is None:
        args.__dict__['output'] = sys.stdout
    else:
        args.__dict__['output'] = open(args.output, 'w')

    return args

if __name__ == '__main__':
    args = parse_args()
    tmpname = os.path.extsep.join(args.bw_in.split(os.path.extsep))[:-1] + '.tmp'
    bw2bed(args.bw_in, parse_chrlen(args.chr_len_file), tmpname, args.output)
