import os
import re
import shlex as sh
import subprocess as sp
from collections import OrderedDict


from config import BOWTIE_EXEC, SAMTOOLS_EXEC


def parse_bowtie_stats(bt_stats):
    # parse this output:
    # <line>*
    # 10263 reads;
    # of these:
    #     10263(100.00 %) were unpaired;
    #     of these:
    #         1055(10.28 %) aligned 0 times
    #         4337(42.26 %) aligned exactly 1 time
    #         4871(47.46 %) aligned > 1 times
    # 89.72 % overall alignment rate
    # <line>*
    N = OrderedDict()
    for i, line in enumerate(bt_stats):
        m = re.match('(\d+) reads.*', line)
        if m is not None:
            N['total'] = m.group(1)
            continue
        m = re.match('(\d)+.*\s+aligned 0 times.*', line)
        if m is not None:
            N['unaligned'] = m.group(1)
            continue
        m = re.match('(\d)+.*\s+aligned exactly 1 time.*', line)
        if m is not None:
            N['unique-align'] = m.group(1)
            continue
        m = re.match('(\d)+.*\s+aligned > 1 times.*', line)
        if m is not None:
            N['multiple-align'] = m.group(1)
            continue
    if len(N) < 4: raise ValueError('\n'.join(bt_stats))
    return N

#
#
# def align(**kwargs):
#     out = sp.PIPE if not kwargs['cnt'] else open(os.devnull)
#     bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (BOWTIE_EXEC, kwargs['n'], kwargs['fastq'], kwargs['index'])),
#                   stdout=out, stderr=sp.PIPE)
#     if not kwargs['cnt']:
#         if kwargs['sam_header_file'] is not None:
#             awkcmd = ''.join(("""awk '{if (substr($1,1,1) == "@" && substr($2,1,2) == "SN")""",
#                               """{print $0 > "%s";} print; }' """)) % kwargs['sam_header_file']
#             bt = sp.Popen(sh.split(awkcmd), stdin=bt.stdout, stdout=sp.PIPE)
#         st = sp.Popen(sh.split('%s view -b -o %s' % (SAMTOOLS_EXEC,kwargs['tmp_bam'])), stdin=bt.stdout)
#         st.wait()
#         if kwargs['unaligned_bam'] is not None:
#             cmd = 'samtools view -f4 -b %s -o %s' % (kwargs['tmp_bam'], kwargs['unaligned_bam'])
#             naligned = sp.Popen(sh.split(cmd), stdout=sp.PIPE)
#             naligned.wait()
#         aligned = sp.Popen(sh.split('%s view -F4 -b %s' % (SAMTOOLS_EXEC, kwargs['tmp_bam'])), stdout=sp.PIPE)
#         sort = sp.Popen(sh.split('%s sort -o %s -T %s' % (SAMTOOLS_EXEC, kwargs['bam'], kwargs['tmp_bam'])),
#                         stdin=aligned.stdout)
#         sort.wait()
#         os.remove(kwargs['tmp_bam'])
#     bt.wait()
#     Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
#     with open(kwargs['align_stats_file'], 'w') as S:
#         S.write('\n'.join('%s %s' % (k,v) for k,v in Ns.items()))


def single_end_count(fastq: str, index: str, n: int=4) -> dict:
    bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (BOWTIE_EXEC, n, fastq, index)),
                  stdout=open(os.devnull), stderr=sp.PIPE)
    bt.wait()
    return parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))


def single_end_align(fastq:str, index:str, bam: str, sam_header:str=None, unaligned_bam:str=None, n:int=4) -> dict:
    tmp_bam = bam + '.tmp'
    bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (BOWTIE_EXEC, n, fastq, index)),
                  stdout=sp.PIPE, stderr=sp.PIPE)
    if sam_header is not None:
            awkcmd = ''.join(("""awk '{if (substr($1,1,1) == "@" && substr($2,1,2) == "SN")""",
                              """{print $0 > "%s";} print; }' """)) % sam_header
            bt = sp.Popen(sh.split(awkcmd), stdin=bt.stdout, stdout=sp.PIPE)
    st = sp.Popen(sh.split('%s view -b -o %s' % (SAMTOOLS_EXEC, tmp_bam)), stdin=bt.stdout)
    st.wait()
    if unaligned_bam is not None:
            cmd = 'samtools view -f4 -b %s -o %s' % (tmp_bam, unaligned_bam)
            naligned = sp.Popen(sh.split(cmd), stdout=sp.PIPE)
            naligned.wait()
    aligned = sp.Popen(sh.split('%s view -F4 -b %s' % (SAMTOOLS_EXEC, tmp_bam)), stdout=sp.PIPE)
    sort = sp.Popen(sh.split('%s sort -o %s -T %s' % (SAMTOOLS_EXEC, bam, tmp_bam)),
                    stdin=aligned.stdout)
    sort.wait()
    os.remove(tmp_bam)
    bt.wait()
    return parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))