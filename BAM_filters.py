
import os
import subprocess as sp
from abc import ABCMeta, abstractmethod
import shlex as sh


class FilterScheme(object):

    def __init__(self, name, filters=None):
        if not filters: filters = []
        self.filters = filters
        self.name = name

    def filter(self, bamin, bamout, sample):
        """
        :param bamin: input bam path
        :param bamout: output bam path
        :return: the number of alignments that passed the filterscheme
        """
        fin = open(bamin, 'rb')
        for f in self.filters:
            fout = f.filter(fin, sample)
            fin = fout
        cnt = sp.Popen(sh.split("""awk 'END {print NR > "/dev/stderr"}'"""),
                       stderr=sp.PIPE, stdout=sp.PIPE, stdin=fin)
        final = sp.Popen(sh.split('samtools sort'), stdin=cnt.stdout, stdout=open(bamout, 'wb'))
        final.wait()
        n = int(''.join(cnt.stderr.read().decode('utf-8')))
        index = sp.Popen(sh.split('samtools index %s' % bamout))
        index.wait()
        return n

    def __str__(self):
        return '\n'.join(str(self.filters))


class BAMFilter(object):
    __metaclass__ = ABCMeta
    name = None
    description = None
    args = None

    def __init__(self, negate, **kwargs):
        self.negate = negate
        self.__dict__.update(kwargs)

    @abstractmethod
    def filter(self, fin, sample):
        """
        filter the data in fin into fout. Input can be assumed to be sorted by genomic position, and output should
        remain sorted similarly.

        :param fin:  filehandle for input
        :return: filehandle for output
        """
        pass

    def __str__(self):
        argstr = ','.join('%s:%s' % (k, str(v)) for k,v in self.__dict__.items())
        return 'negate: %s, %s' % (str(self.negate), argstr)


class DuplicateFilter(BAMFilter):
    name = 'dup'
    description = 'remove artificially amplified reads'
    args = {'kind': (str, 'the method of choice for removal: "start&umi", "start&umi&cigar"')}

    def filter(self, fin, sample):
        hdrpath = os.sep.join([self.bam_dir, sample + '.hdr.sam'])  # ugly, but no choice :(
        inbam = sp.Popen(sh.split('samtools view -S'), stdout=sp.PIPE, stdin=fin)
        ident = '$1'  # umi
        if self.kind == 'start&umi&cigar': ident += '$6' #cigar string
        if self.negate:
            awkstr = ''.join([""" {if (prev!=$4) {delete u; u[%s]="";} """,
                              """ else {if (%s in u) print $0; else u[%s]="";} prev=$4;} """]) % (ident, ident, ident)
        else:
            awkstr = ''.join([""" {if (prev != $4) {for (x in u) {print u[x];} delete u; u[%s]=$0;} """,
                              """ else {u[%s]=$0;} prev=$4;} """,
                              """ END {for (x in u) print u[x];} """]) % (ident, ident)
        awk = sp.Popen(["awk", awkstr], stdin=inbam.stdout, stdout=sp.PIPE)
        addh = sp.Popen(sh.split('cat %s -' % hdrpath), stdin=awk.stdout, stdout=sp.PIPE)
        st = sp.Popen(sh.split('samtools view -b'), stdin=addh.stdout, stdout=sp.PIPE)
        return st.stdout


class AlignmentQualityFilter(BAMFilter):
    name = 'qual'
    description = 'remove reads with low alignment quality (see bowtie/SAM documentation)'
    args = {'q': (int, 'the minimal quality for an alignment to pass the filter')}

    def filter(self, fin, sample):
        if self.negate:
            hdrpath = os.sep.join([self.bam_dir, sample + '.hdr.sam']) # ugly, but no choice :(
            st1 = sp.Popen(sh.split('samtools view -S'), stdout=sp.PIPE, stdin=fin)
            awk = sp.Popen(sh.split("awk '{if ($5 <= %i) {print;}}'" % self.q), stdin=st1.stdout, stdout=sp.PIPE)
            addh = sp.Popen(sh.split('cat %s -' % hdrpath), stdin=awk.stdout, stdout=sp.PIPE)
            st = sp.Popen(sh.split('samtools view -b'), stdin=addh.stdout, stdout=sp.PIPE)
        else:
            st = sp.Popen(sh.split('samtools view -q %i' % self.q), stdin=fin, stdout=sp.PIPE)
        return st.stdout


class PolyAFilter(BAMFilter):
    name = 'polyA'
    description = 'keep only reads that are show evidence of polyadenylation - alignment issues and dA or dT tracks'
    args = {'n': (int, "minimal number of A/T at 3' end to consider a read for a polyA read")}

    def filter(self, fin, sample):
        return fin
        # input = sp.Popen(sh.split('samtools view'), stdout=sp.PIPE, stdin=fin)
        # ident = '$1'  # umi
        # if self.kind == 'start&umi&cigar': ident += '$6'  # cigar string
        # if self.negate:
        #     awkstr = ''.join([""" {if (prev!=$4) {delete u; u[%s]="";} """,
        #                       """ else {if (%s in u) print $0; else u[%s]="";} prev=$4;} """]) % (ident, ident, ident)
        # else:
        #     awkstr = ''.join([""" {if (prev != $4) {for (x in u) {print u[x];} delete u; u[%s]=$0;} """,
        #                       """ else {u[%s]=$0;} prev=$4;} """,
        #                       """ END {for (x in u) print u[x];} """]) % (ident, ident)
        # awk = sp.Popen(['awk', awkstr], stdin=input.stdout, stdout=sp.PIPE)
        # st = sp.Popen(sh.split('samtools view -bq %i'), stdin=awk.stdout, stdout=sp.PIPE)
        # return st.stdout
