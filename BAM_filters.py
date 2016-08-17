
import subprocess as sp
from abc import ABCMeta, abstractmethod
import shlex as sh


class FilterScheme(object):

    def __init__(self, name, filters=None):
        if not filters: filters = []
        self.filters = filters
        self.name = name

    def filter(self, bamin, bamout):
        """
        :param bamin: input bam path
        :param bamout: output bam path
        :return: the number of alignments that passed the filterscheme
        """
        fin = open(bamin, 'rb')
        for f in self.filters:
            fout = f.filter(fin)
            fin = fout
        cnt = sp.Popen(sh.split("""awk 'END {print NR > "/dev/stderr"}'"""),
                       stderr=sp.PIPE, stdout=sp.PIPE, stdin=fin)
        final = sp.Popen(sh.split('samtools sort'), stdin=cnt.stdout, stdout=open(bamout, 'wb'))
        final.wait()
        n = int(''.join(cnt.stderr.read().decode('utf-8')))
        index = sp.Popen(sh.split('samtools index %s' % bamout))
        index.wait()
        return n


class BAMFilter(object):
    __metaclass__ = ABCMeta
    name = None
    description = None
    args = None

    def __init__(self, negate, **kwargs):
        self.negate = negate
        self.__dict__.update(kwargs)

    @abstractmethod
    def filter(self, fin):
        """
        filter the data in (path) bamin into bamout. Sorting and indexing of output is
        performed by wrapping calls, so do not waste time on this.

        :param fin:  filehandle for input
        :return: filehandle for output
        """
        pass


class DuplicateFilter(BAMFilter):
    name = 'dup'
    description = 'remove artificially amplified reads'
    args = {'type': (str, 'the method of choice for removal: "umi&start", "umi&end", "umi&start&end"')}

    def write_awk_script(self):
        pass

    def filter(self, fin):
        return fin


class AlignmentQualityFilter(BAMFilter):
    name = 'qual'
    description = 'remove reads with low alignment quality (see bowtie/SAM documentation)'
    args = {'q': (int, 'the minimal quality for an alignment to pass the filter')}

    def filter(self, fin):
        if self.negate:
            st1 = sp.Popen(sh.split('samtools view'), stdout=sp.PIPE, stdin=fin)
            awk = sp.Popen(sh.split("'{if ($5 <= %i) { print; } }'" % self.q), stdin=st1.stdout, stdout=sp.PIPE)
            st = sp.Popen(sh.split('samtools view -b'), stdout=sp.PIPE)
        else:
            st = sp.Popen(sh.split('samtools view -bq %i'), stdin=fin, stdout=sp.PIPE)
        return st.stdout


class PolyAFilter(BAMFilter):
    name = 'polyA'
    description = 'keep only reads that are show evidence of polyadenylation - alignment issues and dA or dT tracks'
    args = {'n': (int, "minimal number of A/T at 3' end to consider a read for a polyA read")}

    def filter(self, fin):
        return fin
