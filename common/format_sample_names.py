#! /bin/python3
import argparse
import sys


def process(args):
    vars = args.input.readline().strip().split(args.input_delim)[1:]
    seen_s = set([])
    seen_b  = set([])
    alphabet = set('ACGT')
    for i, line in enumerate(args.input):
        line = line.strip().split(args.input_delim)
        s = '_'.join('%s-%s' % vv for vv in zip(vars, line[1:]))
        if s in seen_s: raise ValueError('%s observed twice (at line %i)!' % (s, i))
        seen_s.add(s)
        if line[0] in seen_b: raise ValueError('%s observed twice (at line %i)!' % (line[0], i))
        seen_s.add(line[0])
        if len(set(line[0]) | alphabet) > len(alphabet):
            raise ValueError('first column must be a barcode of %s (error at line %i, bc: '
                             '%s)' % (','.join(alphabet), i, line[0]))
        args.output.write('%s%s%s\n' % (line[0], args.output_delim, s))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input_delim', '-id', type=str, help='delimiter in input', default=',')
    p.add_argument('--output_delim', '-od', type=str, help='delimiter for output', default=',')
    p.add_argument('--input', '-i', type=str, help='input file, default: stdin', default=None)
    p.add_argument('--output', '-o', type=str, help='output file, default: stdout', default=None)
    args = p.parse_args()
    args.__dict__['output'] = open(args.output, 'w') if args.output else sys.stdout
    args.__dict__['input'] = open(args.input, 'rU') if args.input else sys.stdin
    return args


if __name__ == '__main__':
    process(parse_args())