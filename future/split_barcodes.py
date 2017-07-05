import argparse
import subprocess as sp
import os
import shlex as sh
import sys

"""
usage example - split the file pair to barcodes in file "sample_barcodes" with the barcode composed from first 3
bases of R1 and 4,5,6th bases of R2. Perform this in 3 passes to improve performance (remove exact matches quickly, and
then deal with errors), writing to different folders along the way.

$ paste <(zcat fastq/1_S1_R1_001.fastq.gz) <(zcat fastq/1_S1_R2_001.fastq.gz) | paste - - - - | head -100000 |
  python ~/Dropbox/workspace/seqtools/future/split_barcodes.py sample_barcodes.csv -bcf "R1:1,3+R2:4,3" -od hd0 |
  python ~/Dropbox/workspace/seqtools/future/split_barcodes.py sample_barcodes.csv -bcf "R1:1,3+R2:4,3" -d 1 -od hd1 |
  python ~/Dropbox/workspace/seqtools/future/split_barcodes.py sample_barcodes.csv -bcf "R1:1,3+R2:4,3" -d 2 -od hd2 > QC/no-barcode.olf

To merge results, and write to R1/R2 files, or any other modified fastq format (if iterating over samples with variable "s"):
$ cat hd*/$s* | awk -F "\t" '{print $1"\n"$3"\n+\n"$7}' | gzip > fastq/$s.R1.fastq.gz
$ cat hd*/$s* | awk -F "\t" '{print $2"\n"$4"\n+\n"$8}' | gzip > fastq/$s.R2.fastq.gz
"""

NO_BARCODE = 'no-barcode'
FMAP = {        #mapping content to awk index
    'H1': 1,    # header of read 1
    'H2': 2,    # header of read 2
    'R1': 3,    # read 1
    'R2': 4,    # read 2
    'P1': 5,    # plus symbol
    'P2': 6,    # plus symbol
    'Q1': 7,    # quality score of read 1
    'Q2': 8     # quality score of read 2
}


def hamming_ball(w, d, alphabet='ACGTN'):
    """
    A recursive method that returns all words in a hamming ball around w.
    e.g. hamming_ball('aa',1,'ab') -> ['aa','ab','ba']

    :param w: word around which set is generated (center of ball)
    :param d: maximal distance between w and every word in output (radius of ball)
    :param alphabet: list of characters that define the word space
    :return: a set of words in a hamming distance <= d from w in the space defined by alphabet
    """
    ball = [w]
    if d > 0:
        for i in range(len(w)):
            for l in alphabet:
                wl = list(w)
                wl[i] = l
                ball.extend(hamming_ball(''.join(wl), d-1, alphabet))
    return list(set(ball))


def bc_compiler(formatstr):
    bcf = []
    for part in formatstr.split('+'):
        r, params = part.split(':')
        fr, l = params.split(',')
        bcf.append('substr($%i,%i,%i)' % (FMAP[r], int(fr), int(l)))
    return ''.join(bcf)


def compile_awk_command(outdir, cnt_file, bc2sample, bcformat_str): #TODO add barcode parsing from reads
    # awk script will:
    # print every line to a file at outdir/sample.olf according to read barcode
    # print unassigned reads to stdout
    # print barcode counts to cnt_file
    arraydef = ';'.join('a["%s"]="%s.olf"' % (b, s) for b, s in bc2sample.items()) + ';'
    awk_str = ("""'BEGIN {%s} {x=%s; if (x in a) """ % (arraydef, bc_compiler(bcformat_str)),
               """{c[a[x]]++; if (c[a[x]] == 1) printf("") > "%s/"a[x]; """ % outdir,
               """print >> "%s/"a[x];} else {c["%s"]++; print;} }""" % (outdir, NO_BARCODE),
               """END { print "sample","count" > "%s"; for (bc in c) print bc, c[bc] >> "%s" } '""" % (cnt_file, cnt_file))
    return ''.join(awk_str)


def read_barcodes_file(bc_file):
    b2s, snames = {}, set([])
    with open(bc_file) as BCF:
        for line in BCF:
            bc, sname = line.strip().split(',')
            if '_' in sname: raise TypeError('sample name must not contain underscores')
            if bc in b2s: raise ValueError('barcode must be unique (%s is not)' % bc)
            if sname in snames: raise ValueError('sample names must be unique (%s is not)' % sname)
            snames.add(sname)
            b2s[bc] = sname
    return b2s


def split_barcodes(args):
    b2s = read_barcodes_file(args.barcodes)
    e_b2s = {}
    for b, s in b2s.items():
        e_b2s.update({eb: s for eb in hamming_ball(b, args.hamming_distance)})
    b2s = e_b2s

    awkcom = compile_awk_command(args.output_dir, args.count_file_path, b2s, args.bc_format)
    if not os.path.isdir(args.output_dir): os.mkdir(args.output_dir)
    awk = sp.Popen(sh.split('awk -F "\\t" ' + awkcom), stdin=args.input, stdout=args.output)
    awk.wait()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('barcodes', type=str, help='a path to a csv file containing <barcode,sample_name\\n> entries')
    p.add_argument('--input', '-i', default=None,
                   help='From which input is read. Expected as in "FMAP" at top of this file. Default is standrd input.'
                        'to pass to stdin, simply:'
                        '$ paste <(zcat R1.fastq.gz) <(zcat R2.fastq.gz) | paste - - - - | <this script>')
    p.add_argument('--output', '-o', default=None,
                   help='To which main output (orphan reads) is written, in one-liner-fastq format, as in input. '
                        'Default is standard output')
    p.add_argument('--output_dir', '-od', default='.',
                   help='To which output files are written, one per barcode. data is appanded to existing files.')
    p.add_argument('--bc_format', '-bcf', type=str, default='R2:1,7',
                   help='The bases that compile the expected barcode, when concatenated, e.g. R1:1,7+R2:20,3 means that'
                        'the first 7 bases of R1 are concatented to 3 bases (20/21/22) of R2 to create a 10bp sample'
                        ' barcode')
    p.add_argument('--count_file_name', '-cf', default='barcode_counts.csv', type=str,
                   help='To which barcode counts are saved, within output folder.')
    p.add_argument('--hamming_distance', '-d', default=0, type=int,
                   help='Reads will be assigned to samples if their barcode is within this hamming distance from '
                        'respective barcode. If overlap arises, assignment behavior is not guaranteed.')

    args = p.parse_args()
    args.__dict__['output'] = sys.stdout if args.output is None else open(args.output, 'w')
    args.__dict__['input'] = sys.stdin if args.input is None else open(args.input)
    args.__dict__['count_file_path'] = args.output_dir + os.path.sep + args.count_file_name
    return args


if __name__ == '__main__':
    split_barcodes(parse_args())