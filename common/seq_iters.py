import os, gzip, sys

sys.path.append(os.path.abspath('/cs/wetlab/Alon/Dropbox/workspace/'))
sys.path.append(os.path.abspath('/Users/user/Dropbox/workspace/'))

from utils.general import Object
def canonic_path(fname): return os.path.abspath(os.path.expanduser(fname))

DNA2BISULF = dict(N='N',G='G',C='T',T='T',A='A')
REVCOMP = lambda s: ''.join(DNA[b] for b in reversed(s))
TRANS = lambda s, fr: ''.join([GCODE[s[i:i+3]] for i in range(fr,len(s)-3,3)])
DNA = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
SAM_FIELDS = Object(id=0, flag=1, ref=2, start=3, qual=6, CIGAR=5)

STOP = '$'

GCODE = {'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
         'TTA' : 'L', 'TTG' : 'L', 'CTT' : 'L', 'CTC' : 'L',
         'CTA' : 'L', 'CTG' : 'L', 'CGT' : 'R', 'CGC' : 'R',
         'CGA' : 'R', 'CGG' : 'R', 'AGA' : 'R', 'AGG' : 'R',
         'AAA' : 'K', 'AAG' : 'K', 'AAT' : 'N', 'AAC' : 'N',
         'ATG' : 'M', 'GAT' : 'D', 'GAC' : 'D', 'TTT' : 'F',
         'TTC' : 'F', 'TGT' : 'C', 'TGC' : 'C', 'CCT' : 'P',
         'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P', 'CAA' : 'Q', 
         'CAG' : 'Q', 'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S',
         'TCG' : 'S', 'AGT' : 'S', 'AGC' : 'S', 'GAA' : 'E',
         'GAG' : 'E', 'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T',
         'ACG' : 'T', 'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G',
         'GGG' : 'G', 'TGG' : 'W', 'CAT' : 'H', 'CAC' : 'H',
         'TAT' : 'Y', 'TAC' : 'Y', 'ATT' : 'I', 'ATC' : 'I',
         'ATA' : 'I', 'GTT' : 'V', 'GTC' : 'V', 'GTG' : 'V',
         'GTA' : 'V', 'TAA' : STOP, 'TGA' : STOP, 'TAG' : STOP}


def revcomp(seq, aslist=False):
    seq = [DNA[x] for x in seq[::-1]]
    if not aslist: seq = ''.join(seq)
    return seq

def hamming_ball(seq, radius, alphabet='CGTAN'):
    ball = [seq]
    if radius > 0:
        for i in range(len(seq)):
            for l in alphabet:
                seql = list(seq)
                seql[i] = l
                ball.extend(hamming_ball(''.join(seql), radius-1, alphabet))
    return ball


def fasta_iter(fname):
    '''A generator over a fasta file.
    At each generation, returns a (header, sequence) pair.
    '''
    fhandle = open(canonic_path(fname), 'r')
    line = fhandle.readline()
    while line:
        header = line.rstrip()
        seq = []
        line = fhandle.readline()
        while line and line[0] != '>': #collecting sequence
            seq.append(line.rstrip())
            line = fhandle.readline()
        seq = ''.join(seq)
        yield (header,seq)


def rev_genome(fin,fout):
    fout = open(fout,'wb')
    for hdr,seq in fasta_iter(fin):
        fout.write(hdr+'\n')
        fout.write(REVCOMP(seq)+'\n')
        

SAM_REFSTRAND = 0x0010
SAM_SUBSEQUENT = 0x0100

def SAM_iter(samtext):
    '''iterate over seqs in SAM format, for each seq returns a list of alignemnts.
    An alignment entry consists of:
        
    '''
    stiter = (line for line in samtext)
    i, prevlines= 0, []
    line = stiter.next().rstrip().split('\t')
    while line[0]:
        while line[0] and int(line[0]) == i:
            prevlines.append(line)
            line = stiter.next().rstrip().split('\t')
        yield [(l[2], int(l[3]), int(l[1])&SAM_REFSTRAND==0, l[5]) for l in prevlines]
        i, prevlines = i+1, [line]


def multiple_alignment_iter(samtext):
    '''iterate over seqs in SAM format, for each seq returns a list of alignemnts.
    An alignment entry consists of:
    '''
    alignbatch = None
    for line in samtext:
        if not line: continue
        l = line.rstrip().split('\t')
        flag = int(l[SAM_FIELDS.flag])
        if not (flag & SAM_SUBSEQUENT): #new read, yield previous and new batch
            if alignbatch is not None: yield alignbatch
            alignbatch = []
        if l[2] != '*':
            alignbatch.append(Object(chr=l[2],pos=int(l[3]),strand=int(l[1])&SAM_REFSTRAND==0,cigar=l[5]))
    yield alignbatch


def fastq_iter(fnames=None):
    '''Read generator from a collection of fastq files.
    if files and with gz they are treated as gzipped.

    format example:
    @SBS123:140:C23D0ACXX:1:1101:1168:2177 1:N:0:GAAGAAGT #header
    AGGGAAGTCGGCAAAATAGATCCGTAACTTCGGGATAAGGATTGGCTCN     #sequence
    +
    BBBFFFFFFFFFFIIFFIIBF7BBBFFFIIFFFB0BBFFFFF<FFFIF#     #quality
    '''
    for fname in fnames:
        if fname is None:
            while True: yield ('dummy', 'NNNNN', 'XXXXX')
        zipped = fname.endswith('.gz')
        if zipped:
            fhandle = gzip.GzipFile(canonic_path(fname), 'r')
        else:
            fhandle = open(canonic_path(fname), 'r')
        line = fhandle.readline()
        while line:
            header = line.rstrip()
            seq = fhandle.readline().rstrip()
            fhandle.readline()  # +
            qual = fhandle.readline().rstrip()
            yield (header, seq, qual)
            line = fhandle.readline()