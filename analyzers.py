
import os
import subprocess as sp
from abc import ABCMeta, abstractmethod
import shlex as sh
import sys

from config import *
from config import SCER_GENOME_LENGTH_PATH as SGLP
from utils import *


def collect_analyzers():
    analyzers = {}
    analyze_module = sys.modules[__name__]
    for x in dir(analyze_module):
        try:
            if issubclass(analyze_module.__dict__[x], BAMAnalyzer):
                aname = analyze_module.__dict__[x].name
                if aname is not None:
                    analyzers[aname] = analyze_module.__dict__[x]
        except TypeError: pass
    return analyzers


class BAMAnalyzer(object):
    __metaclass__ = ABCMeta
    name = None
    description = None
    args = None

    def __init__(self, tmp_dir, **kwargs):
        self.tmp_dir = tmp_dir
        self.__dict__.update(kwargs)

    @abstractmethod
    def analyze(self, name, bamin, outpath):
        """

        :param bamin: input bam (sorted, indexed) file
        :param outpath: output path, does not include the file names
        :return: a 2-tuple:
            0: a dictionary containing all generated files, and their types (name : type) that will be used to report
               progression and are passed to the aggregate function for easy aggregation access
            1:  count information regarding this sample (QC related). e.g. total # reads used by analysis.
                count format is a dictionary with (property : #) entries. These are added to the general run statistics.
        """
        pass

    def __str__(self):
        return self.name


class Track3p(BAMAnalyzer):
    name = "3pT"
    description = "make strand-specific 3' count tracks"
    args = {}
    suffixes = {'.w.3p.bw': ('bigwig', "watson strand 3' counts"),
                '.c.3p.bw': ('bigwig', "crick strand 3' counts")}
    out_folder = "3p_tracks"

    def analyze(self, files):
        tmp = self.tmp_dir + os.sep + os.path.split(files['.w.3p.bw'])[1].split(os.extsep)[0] + '.3p.bed.tmp'
        def handle_strand(char, fout):
            bedcmd = "bedtools genomecov -3 -ibam %s -g %s  -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam_in'], SGLP, STRANDS[char])), stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(tmp,'w'))
            sbed.wait()
            bw = sp.Popen([BG2W_EXEC, tmp, SGLP, fout])
            bw.wait()
            os.remove(tmp)
        handle_strand('w', files['.w.3p.bw'])
        handle_strand('c', files['.c.3p.bw'])
        return {}


class TrackCoverage(BAMAnalyzer):
    name = "covT"
    description = "standard coverage tracks"
    args = {}
    suffixes = {'.w.bw': ('bigwig', "watson strand coverage"),
                '.c.bw': ('bigwig', "crick strand coverage")}
    out_folder = "coverage_tracks"

    def analyze(self, files):
        tmp = self.tmp_dir + os.sep + os.path.split(files['.w.bw'])[1].split(os.extsep)[0] + '.cov.bed.tmp'
        def handle_strand(char, fout):
            bedcmd = "bedtools genomecov -ibam %s -g %s -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam_in'], SGLP, STRANDS[char])), stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(tmp, 'w'))
            sbed.wait()
            bw = sp.Popen([BG2W_EXEC, tmp, SGLP, fout])
            bw.wait()
            os.remove(tmp)
        handle_strand('w', files['.w.bw'])
        handle_strand('c', files['.c.bw'])
        return {}


class Track5p(BAMAnalyzer):
    name = "5pT"
    description = "make strand-specific 5' count tracks"
    args = {}
    suffixes = {'.w.5p.bw': ('bigwig', "watson strand 5' counts"),
                '.c.5p.bw': ('bigwig', "crick strand 5' counts")}
    out_folder = "5p_tracks"

    def analyze(self, files):
        tmp = self.tmp_dir + os.sep + os.path.split(files['.w.5p.bw'])[1].split(os.extsep)[0] + '.5p.bed.tmp'
        def handle_strand(char, fout):
            bedcmd = "bedtools genomecov -5 -ibam %s -g %s  -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam_in'], SGLP, STRANDS[char])), stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(tmp,'w'))
            sbed.wait()
            bw = sp.Popen([BG2W_EXEC, tmp, SGLP, fout])
            bw.wait()
            os.remove(tmp)
        handle_strand('w', files['.w.5p.bw'])
        handle_strand('c', files['.c.5p.bw'])
        return {}


class TrackCoverage(BAMAnalyzer):
    name = "covT"
    description = "standard coverage tracks"
    args = {}
    suffixes = {'.w.bw': ('bigwig', "watson strand coverage"),
                '.c.bw': ('bigwig', "crick strand coverage")}
    out_folder = "coverage_tracks"

    def analyze(self, files):
        tmp = self.tmp_dir + os.sep + os.path.split(files['.w.bw'])[1].split(os.extsep)[0] + '.cov.bed.tmp'
        def handle_strand(char, fout):
            bedcmd = "bedtools genomecov -ibam %s -g %s -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam_in'], SGLP, STRANDS[char])), stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(tmp, 'w'))
            sbed.wait()
            bw = sp.Popen([BG2W_EXEC, tmp, SGLP, fout])
            bw.wait()
            os.remove(tmp)
        handle_strand('w', files['.w.bw'])
        handle_strand('c', files['.c.bw'])
        return {}
