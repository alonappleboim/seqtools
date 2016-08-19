
import os
import subprocess as sp
from abc import ABCMeta, abstractmethod
import shlex as sh
import sys


def collect_filters():
    filters = {}
    filter_module = sys.modules[__name__]
    for x in dir(filter_module):
        try:
            if issubclass(filter_module.__dict__[x], SAMFilter):
                fname = filter_module.__dict__[x].name
                if fname is not None:
                    filters[fname] = filter_module.__dict__[x]
        except TypeError: pass
    return filters


class FilterScheme(object):

    def __init__(self, name, filters=None):
        if not filters: filters = []
        self.filters = filters
        self.name = name

    def filter(self, bamin, bamout, sam_hdr):
        """
        :param bamin: input bam path
        :param bamout: output bam path
        :param sam_hdr: header associated with bam input
        :return: the number of alignments that passed the filter scheme
        """
        fin = sp.Popen(sh.split('samtools view'), stdout=sp.PIPE, stdin=open(bamin, 'rb')).stdout
        for f in self.filters:
            fout = f.filter(fin)
            fin = fout
        addh = sp.Popen(sh.split('cat %s -' % sam_hdr), stdin=fin, stdout=sp.PIPE)
        bam = sp.Popen(sh.split('samtools view -b'), stdout=sp.PIPE, stdin=addh.stdout)
        final = sp.Popen(sh.split('samtools sort'), stdin=bam.stdout, stdout=open(bamout, 'wb'))
        final.wait()
        index = sp.Popen(sh.split('samtools index %s' % bamout))
        index.wait()
        cnt = sp.Popen(sh.split("samtools idxstats %s" % bamout), stdout=sp.PIPE)
        out = ''.join(cnt.communicate()[0].decode('utf-8'))
        n = 0
        for line in out.split('\n'):
            if not line or line.startswith('*'): continue
            n += int(line.split('\t')[2])
        return n

    def __str__(self):
        return '\n'.join(['%s| %s' % (x.name, str(x)) for x in self.filters])


class SAMFilter(object):
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
        filter the data in fin into fout. Input can be assumed to be sorted by genomic position, and output should
        remain sorted similarly.

        :param fin:  filehandle for SAM input without a header
        :return: filehandle for SAM output without a header
        """
        pass

    def __str__(self):
        return ','.join('%s:%s' % (k, str(v)) for k,v in self.__dict__.items())


class DuplicateFilter(SAMFilter):
    name = 'dup'
    description = 'remove artificially amplified reads'
    args = {'kind': (str, 'start&umi&cigar', 'the method of choice for removal: "start&umi", "start&umi&cigar"')}

    def filter(self, fin):
        ident = '$1'  # umi
        if self.kind == 'start&umi&cigar': ident += '$6' #cigar string
        if self.negate:
            awkstr = ''.join([""" {if (prev!=$4) {delete u; u[%s]="";} """,
                              """ else {if (%s in u) print $0; else u[%s]="";} prev=$4;} """]) % (ident, ident, ident)
        else:
            awkstr = ''.join([""" {if (prev != $4) {for (x in u) {print u[x];} delete u; u[%s]=$0;} """,
                              """ else {u[%s]=$0;} prev=$4;} """,
                              """ END {for (x in u) print u[x];} """]) % (ident, ident)
        f = sp.Popen(["awk", awkstr], stdin=fin, stdout=sp.PIPE)
        return f.stdout


class AlignmentQualityFilter(SAMFilter):
    name = 'qual'
    description = 'remove reads with low alignment quality (see bowtie/SAM documentation)'
    args = {'q': (int, 5, 'the minimal quality for an alignment to pass the filter')}

    def filter(self, fin):
        comparator = '<=' if self.negate else '>'
        f = sp.Popen(sh.split("awk '{if ($5 %s %i) {print;}}'" % (comparator, self.q)),
                     stdin=fin, stdout=sp.PIPE)
        return f.stdout


class PolyAFilter(SAMFilter):
    name = 'polyA'
    description = 'keep only reads that are show evidence of polyadenylation - alignment issues and dA or dT tracks'
    args = {'n': (int, 5, "minimal number of A/T at 3' end to consider a read for a polyA read")}

    def filter(self, fin):
        return fin
        input = sp.Popen(sh.split('samtools view'), stdout=sp.PIPE, stdin=fin)
        ident = '$1'  # umi
        if self.kind == 'start&umi&cigar': ident += '$6'  # cigar string
        if self.negate:
            awkstr = ''.join([""" {if (prev!=$4) {delete u; u[%s]="";} """,
                              """ else {if (%s in u) print $0; else u[%s]="";} prev=$4;} """]) % (ident, ident, ident)
        else:
            awkstr = ''.join([""" {if (prev != $4) {for (x in u) {print u[x];} delete u; u[%s]=$0;} """,
                              """ else {u[%s]=$0;} prev=$4;} """,
                              """ END {for (x in u) print u[x];} """]) % (ident, ident)
        awk = sp.Popen(['awk', awkstr], stdin=input.stdout, stdout=sp.PIPE)
        st = sp.Popen(sh.split('samtools view -bq %i'), stdin=awk.stdout, stdout=sp.PIPE)
        return st.stdout


def scheme_from_parse_tree(name, ptree):
    filters = collect_filters()
    fs = []
    for f, attr in ptree.items():
        if f not in filters:
            raise ValueError('No filter %s' % f)
        F = filters[f]
        args = {}
        for aname, aval in attr['args'].items():
            args[aname] = F.args[aname][0](aval)
        assert attr['neg'] in ['','+','-']
        for a, props in F.args.items():
            if a not in args: args[a] = props[1]  # default
        fs.append(F(negate=attr['neg'] == '-', **args))
    return FilterScheme(name, fs)


def build_filter_schemes(filter_scheme_string):
    pt = parse_filterscheme_list(filter_scheme_string)
    schemes = {}
    for fs in pt.keys():
        schemes[fs] = scheme_from_parse_tree(fs, pt[fs])
    return schemes


def parse_filterscheme_list(fstr):
    """
        format is:
        [<filter_scheme>:[<filter_name>([<argname1=argval1>,]+)[+|-])]+;]*
        examples:
            non-unique:qual(q=10)-,dup()+;unique:dup(),qual(q=10)+;unique-polya:dup()+,qual(q=1)+,polya(nA=5)+

        :return: a dictionary of instantiated filter schemes
        """
    ptree = {}
    for fs in fstr.split(';'):
        [name, pt] = parse_filterscheme_string(fs.strip())
        ptree[name] = pt
    return ptree


def parse_filterscheme_string(fstr):
    name, rest = fstr.strip().split(':')
    ptree = {}
    while True:
        if '(' not in rest: break
        fname, rest = rest.split('(', 1)
        ptree[fname.strip()] = {}
        argstr, rest = rest.split(')', 1)
        if ',' not in rest: neg = rest
        else: neg, rest = rest.split(',', 1)
        ptree[fname]['neg'] = neg.strip()
        args = {}
        for arg in argstr.strip().split(','):
            if arg:
                (aname, aval) = arg.split('=')
                args[aname.strip()] = aval.strip()
        ptree[fname]['args'] = args
    return name, ptree


if __name__ == '__main__':
    print(parse_filterscheme_list('non-unique:qual(q=10)-,dup()+;unique:dup(),qual(q=10)+;unique-polya:dup()+,qual(q=1)+,polyA(n=5)+'))
