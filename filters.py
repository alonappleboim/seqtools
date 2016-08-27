
import os
import copy
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


class FilterPipe(object):
    """
    Collects the "And" result of multiple filters
    """

    def __init__(self, name, filters=None):
        if not filters: filters = []
        self.filters = filters
        self.name = name

    def filter(self, samin, bamout, sam_hdr, bamnfilter=None):
        """
        :param bamin: input bam path
        :param bamout: output bam path
        :param sam_hdr: header associated with bam input
        :param bamnfilter: if given, filtered reads are written to this file
        :return: the number of alignments that passed the filter scheme
        """
        fin = sp.Popen(sh.split('samtools view'), stdout=sp.PIPE, stdin=open(samin, 'rb')).stdout
        tmp_files = []
        for i, f in enumerate(self.filters):
            if bamnfilter is not None:  # split output and pass (also) to negated filter
                tmprep = bamnfilter + '.tmp.' + str(i)
                tmp_files.append(tmprep)
                fin = sp.Popen(['tee', tmprep], stdin=fin, stdout=sp.PIPE).stdout
            fout = f.filter(fin)
            fin = fout
        addh = sp.Popen(sh.split('cat %s -' % sam_hdr), stdin=fin, stdout=sp.PIPE)
        bam = sp.Popen(sh.split('samtools view -b'), stdout=sp.PIPE, stdin=addh.stdout)
        final = sp.Popen(sh.split('samtools sort -T %s' % bamout), stdin=bam.stdout, stdout=open(bamout, 'wb'))
        final.wait()
        index = sp.Popen(sh.split('samtools index %s' % bamout))
        index.wait()
        # TODO: collect duplicate count distribution
        if bamnfilter is not None: #TODO: collect filtration statistics
            with open(bamnfilter+ '.tmp', 'wb') as fout: #TODO: write each filtered stream to designated folder
                for i, (f, tmpf) in enumerate(zip(self.filters, tmp_files)):
                    nf = copy.copy(f)
                    nf.negate = ~nf.negate
                    fout.write(nf.filter(open(tmpf)).read())
                    os.remove(tmpf)
            addh = sp.Popen(['cat', sam_hdr, bamnfilter + '.tmp'], stdout=sp.PIPE)
            sp.Popen(sh.split('samtools view -b'), stdout=open(bamnfilter, 'wb'), stdin=addh.stdout).wait()
            os.remove(bamnfilter + '.tmp')

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
    def filter(self, fin, filtered=None):
        """
        filter the data in fin into fout. Input can be assumed to be sorted by genomic position, and output should
        remain sorted similarly.

        :param fin:  filehandle for (headerless-)SAM input
        :param filtered: a filepath to which filtered entries are written
        :return: filehandle for (headerless-)SAM input
        """
        pass

    def __str__(self):
        return ','.join('%s:%s' % (k, str(v)) for k,v in self.__dict__.items())


class SingleEndFilter(SAMFilter):
    __metaclass__ = ABCMeta


class DuplicateFilter(SingleEndFilter):
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


class AlignmentQualityFilter(SingleEndFilter):
    name = 'qual'
    description = 'remove reads with low alignment quality (see bowtie/SAM documentation)'
    args = {'qmin': (int, 5, 'minimal quality, exclusive'),
            'qmax': (int, 255, 'maximal quality, exclusive')}

    def filter(self, fin):
        awkcmd = "{if ($5 < %i && $5 > %i) {print;}}" % (self.qmax, self.qmin)
        if self.negate:
            awkcmd = "{if ($5 >= %i || $5 <= %i) {print;}}" % (self.qmax, self.qmin)
        f = sp.Popen(['awk', awkcmd], stdin=fin, stdout=sp.PIPE)
        return f.stdout


class PolyAFilter(SingleEndFilter):
    name = 'polyA'
    description = 'keep only reads that show evidence of polyadenylation - alignment issues at edges and dA or dT tracks'
    args = {'n': (int, 6, "minimal number of A/T at 3' end to consider a read for a polyA read"),
            'p': (float, .8, "fraction of extermal unmatched bases that must be A/T")}

    def filter(self, fin):
        awkcmd = """
        {
        c = 0;
        p = 0;
        patsplit($6, T, /[MIDNSHPX]/, N)
        split($10, seq, "");
        SL = length(seq);
        NL = length(N);
        if (and($2,0x16) == 16) {
            if (T[1] == "S" && N[0] > %i) {
                for(i = 1; i <= N[0]; i++) if(seq[i] == "T") c++;
                if (c > (N[0]*%f)) {p = 1;}
            }
        }
        else {
            if (T[NL] == "S" && N[NL-1] > %i) {
                for (i = SL; i > SL - N[NL-1]; i--) if (seq[i] == "A") c++;
                if (c > (N[NL-1]*%f)) {p = 1;}
            }
        }
        if (p == %i) { print;}
        }""".replace("\n","") % (self.n, self.p, self.n, self.p, 0 if self.negate else 1)
        f = sp.Popen(["awk", awkcmd], stdin=fin, stdout=sp.PIPE)
        return f.stdout


class StrandFilter(SingleEndFilter):
    name = 'strand'
    description = 'keep only reads that are in specified strand'
    args = {'s': (str, 'w', "'w' - watson, 'c' - crick")}

    def filter(self, fin):
        if self.s == 'w':
            sym = '==' if not self.negate else '!='
        else:
            sym = '!=' if not self.negate else '=='
        awkcmd = """ { if (and($2,0x16) %s 16) {print;} }""" % sym
        f = sp.Popen(["awk", awkcmd], stdin=fin, stdout=sp.PIPE)
        return f.stdout


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
    return FilterPipe(name, fs)


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
        name, rest = fs.strip().split(':')
        pt = parse_filterscheme_string(rest.strip())
        ptree[name] = pt
    return ptree


def parse_filterscheme_string(rest):
    ptree = {}
    while True:
        if '(' not in rest: break
        fname, rest = rest.split('(', 1)
        ptree[fname.strip()] = {}
        argstr, rest = rest.split(')', 1)
        if ',' not in rest:
            neg = rest
        else:
            neg, rest = rest.split(',', 1)
        ptree[fname]['neg'] = neg.strip()
        args = {}
        for arg in argstr.strip().split(','):
            if arg:
                (aname, aval) = arg.split('=')
                args[aname.strip()] = aval.strip()
        ptree[fname]['args'] = args
    return ptree


if __name__ == '__main__':
    print(parse_filterscheme_list('non-unique:qual(q=10)-,dup()+;unique:dup(),qual(q=10)+;unique-polya:dup()+,qual(q=1)+,polyA(n=5)+'))
