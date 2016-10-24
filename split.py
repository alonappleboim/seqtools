"""
This script performs the preprocessing steps and starts a dedicated process for every sample.
 """


from collections import Counter
import getpass
import csv
import re
from secure_smtp import ThreadedTlsSMTPHandler
import logging as lg
import argparse
import os
import sys
import datetime
import multiprocessing as mp
import subprocess as sp
import shlex as sh
import shutil

from config import *

if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()

from work import *
from utils import *
from filters import *


class Aspect(object):

    def __init__(self, name, shortname, vals, units):
        self.name = name
        self.shortname = shortname
        self.vals = vals
        self.units = units


def cartesian_prod(varvals):
    prod = []
    for var, vals in varvals.items():
        ldim = cartesian({v: vs for v, vs in varvals.items() if v != var})
        for vl in vals:
            prod.append([(var,vl)] + ldim)
    return prod


class Splitter(object):
    """
    Applies a cartesian product of filters to each sample
    [f1, f2] X [f3, f4] = {(f1&f3),(f1&f4),(f2&f3), (f2&f4)}
    """


def split(self, infiles, filter_partition, output_folder):
    """
    :param infiles: prefix -> (bam file, sam_hdr)
    :param filter_partition: a list of filter list that will be multiplied to produce split files:
        {grp1:[(x,filterX),(y,filterY)], grp2:[(a,filterA),(b,filterB)]} =>
        prefix_grp1-x_grp2-a.bam % passed filterX and filterA
        prefix_grp1-y_grp2-a.bam % ...
        prefix_grp1-x_grp2-b.bam
        prefix_grp1-y_grp2-b.bam
    """
    fschemes = {}
    for path in cartesian_prod(filter_partition):
        name = '_'.join('%s-%s' % (g, name) for (g, (name, f)) in path)
        fs = FilterPipe(name, [f for _, (_, f) in path])
        for prefix, (bamin, sam_hdr) in infiles.items():
            bamout = output_folder + os.sep + prefix + '_' + name + '.bam'
            fs.filter(bamin, bamout,sam_hdr)