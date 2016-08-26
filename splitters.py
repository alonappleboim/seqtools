
import os
import subprocess as sp
from abc import ABCMeta, abstractmethod
import shlex as sh
import sys


def collect_splitters():
    splitters = {}
    splitter_module = sys.modules[__name__]
    for x in dir(splitter_module):
        try:
            if issubclass(splitter_module.__dict__[x], SAMFilter):
                fname = filter_module.__dict__[x].name
                if fname is not None:
                    filters[fname] = filter_module.__dict__[x]
        except TypeError: pass
    return splitters


class SplitScheme(object):
    """
    a pipeline of splitters
    """
    def __init__(self, splitters):
        self.splitters = splitters

    def split(self, split_arg_map):


class BAMSplitter():
    __metaclass__ = ABCMeta
    name = None
    aspect_name = None
    description = None
    args = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @abstractmethod
    def split(self, fin, fout, val):
        """
        split the data in fhandles['bamin'] into the files defined by the splitter initialization arguments.

        :param fhandles: a dictionary of file paths
        :return: splitting statistics
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


# class PolyAFilter(SAMFilter):
#     name = 'polyA'
#     description = 'keep only reads that show evidence of polyadenylation - alignment issues at edges and dA or dT tracks'
#     args = {'n': (int, 6, "minimal number of A/T at 3' end to consider a read for a polyA read"),
#             'p': (float, .8, "fraction of extermal unmatched bases that must be A/T")}
#
#     def filter(self, fin):
#         awkcmd = """
#         {
#         c = 0;
#         patsplit($6, T, /[MIDNSHPX]/, N)
#         split($10, seq, "");
#         SL = length(seq);
#         NL = length(N);
#         if (and($2,0x16) == 16) {
#             if (T[1] == "S" && N[0] > %i) {
#                 for(i = 1; i <= N[0]; i++) if(seq[i] == "T") c++;
#                 if (c > (N[0]*%f)) {print;}
#             }
#         }
#         else {
#             if (T[NL] == "S" && N[NL-1] > %i) {
#                 for (i = SL; i > SL - N[NL-1]; i--) if (seq[i] == "A") c++;
#                 if (c > (N[NL-1]*%f)) {print;}
#             }
#         }
#         }""".replace("\n","") % (self.n, self.p, self.n, self.p)
#         f = sp.Popen(["awk", awkcmd], stdin=fin, stdout=sp.PIPE)
#         return f.stdout


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
