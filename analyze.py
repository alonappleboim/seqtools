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


from workers import *
from analyzers import *
from utils import *
from filters import *


class SampleFeature(object):

    def __init__(self, caption, units=None, vals=None, type='str'):
        self.caption = caption
        self.units = units
        self.vals = set(vals)
        self.type = type


class Analyzer(object):

    def setup_log(self):
        logger = lg.getLogger()
        logfile = get_logfile()
        logger.setLevel(lg.DEBUG)
        fh = lg.FileHandler(logfile)  # copy log to the log folder when done.
        fh.formatter = lg.Formatter('%(levelname)s\t%(asctime)s\t%(message)s')
        ch = lg.StreamHandler()
        ch.formatter = lg.Formatter('%(asctime)-15s\t%(message)s')
        fh.setLevel(lg.DEBUG)
        ch.setLevel(lg.INFO)
        if self.debug is not None: ch.setLevel(lg.DEBUG)
        if self.user_emails is not None:
            mailh = ThreadedTlsSMTPHandler(mailhost=('smtp.gmail.com', 587),
                                           fromaddr='transeq.pipeline@google.com',
                                           toaddrs=self.user_emails.split(','),
                                           credentials=('transeq.pipeline', 'transeq1234'),
                                           subject='TRANSEQ pipeline message')
            mailh.setLevel(lg.CRITICAL)
            logger.addHandler(mailh)
        logger.addHandler(fh)
        logger.addHandler(ch)
        self.logger = logger
        self.logfile = logfile

    def __init__(self, argobj, cmdline, samples=None, features=None):
        self.__dict__.update(argobj.__dict__)

        self.setup_log()
        self.logger.log(lg.INFO, 'commandline: %s' % cmdline)
        self.comq = mp.Queue()
        if self.debug:
            self.logger.log(lg.INFO, '=================== DEBUG MODE (%s) ===================' % self.debug)
        self.check_third_party()

        if samples is None or features is None:
            self.files_and_features()
        else:
            self.features = features
            self.samples = samples
        self.generate_dir_tree()

        self.w_manager = WorkManager(self.comq, self.exec_on=='slurm')

    def execute(self):
        self.split_files()
        self.aggregate_reads()
        self.generate_hubs()
        self.generate_annotation_arrays()
        self.export_outputs()
        self.aftermath()

    def split_files(self):
        self.logger.log(lg.INFO, 'Splitting BAM files...')
        token_map = {}
        for files, sample in self.sample_files.items():
            for chain in self.splitgrid:
                subsample = Sample()
                for g, (fname, f) in chain:
                    subsample.feats[g] = fname
                fs = FilterPipe('', [f for _, (_, f) in chain])
                bamout = self.bam_dir + os.sep + subsample.name + '.bam'
                subsample.bam = bamout
                kwargs = dict(bamin=files['bam'], bamout=bamout, sam_hdr=files['sam_hdr'])
                token_map[self.w_manager.run(fs.filter, kwargs)] = subsample
        samples = []
        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            samples.append(sample)
            if e is not None:
                msg = 'error while splitting sample %s:\n%s' % (sample.name, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                stat = 'total' % f
                if stat not in self.stat_order: self.stat_order.append(stat)
                self.logger.log(lg.DEBUG, 'subsample %s was split (%i reads)' % (sample.name, info))
        self.samples = samples
        self.logger.log(lg.INFO, 'Splitting completed.')
        self.mark()

    def aggregate_reads(self):
        self.logger.log(lg.INFO, 'Aggregating...')
        token_map = {}
        for sample in self.samples():
            sample.bw = self.bw_dir + os.sep + sample.name + '.bw'
            args = dict(bamin=sample.bam, bwout=sample.bw)
            token_map[self.w_manager.run(self.aggregator.aggregate, kwargs=args)] = sample
        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error while aggregating %s with %s:\n %s' % (sample, self.aggregator.name, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                self.logger.log(lg.DEBUG, 'aggregated %s with %s.' % (sample, self.aggregator.name))

    def generate_hubs(self):
        self.logger.log(lg.INFO, 'Generating hubs...')
        token_map = {}
        for hub, spec in self.hubs.items():
            args = dict(hubspec=spec, hub_dir=self.hub_dir, samples=self.samples)
            token_map[self.w_manager.run(Analyzer.generate_hub, kwargs=args)] = hub
        while token_map:
            tok, e, info = self.comq.get()
            hub = token_map.pop(tok)
            if e is not None:
                msg = 'error while generating hub %s:\n %s' % (hub, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                self.logger.log(lg.DEBUG, 'generated hub %s: %s.' % (hub, self.hubs[hub].url))
        self.logger.log(lg.INFO, 'Done, see all hubs in %s.' % self.hub_dir)

    def generate_annotation_arrays(self):
        self.logger.log(lg.INFO, 'Generating annotation arrays...')
        token_map = {}
        data_files = {}
        for annot, spec in self.annots.items():
            for sample in self.samples:
                tmp_path = self.tmp_dir + os.sep + annot + '_' +  sample.name + '.annot.cnt'
                sample.tmp_annot_data[annot] = tmp_path
                args = dict(annot_spec=spec, tmp_path=tmp_path, sample=sample)
                token_map[self.w_manager.run(Analyzer.generate_hub, kwargs=args)] = (sample, annot)
        while token_map:
            tok, e, info = self.comq.get()
            sample, annot = token_map.pop(tok)
            if e is not None:
                if annot not in data_files: data_files[annot] = {}
                msg = 'error while extracting annotation data for %s from sample %s:\n %s' % (annot, sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                self.logger.log(lg.DEBUG, 'collected annotation data for %s from sample %s.' % (annot, sample))

        for a, spec in self.annots.items():
            self.output_annot(a)
            self.logger.log(lg.INFO, 'Collected annotation data (%s) into file %s' % (spec.summary, spec.outpath))
        self.logger.log(lg.INFO, 'Collected all data in folder %s.' % self.data_dir)

    def copy_log(self):
        shutil.copy(self.logfile, self.output_dir + os.sep + 'full.log')

    def check_third_party(self):
        if self.exec_on == 'slurm':
            try:
                p = sp.Popen(sh.split('srun "slurm, are you there?"'), stdout=sp.PIPE, stderr=sp.PIPE)
                p.communicate()
                self.logger.log(lg.INFO, "slurm check.. OK")
            except OSError as e:
                self.logger.log(lg.CRITICAL, "This is not a slurm cluster, execute with flag -eo=local")
                raise e

        for ex in [k for k in self.__dict__.keys() if k.endswith('exec')]:
            try:
                p = sp.Popen(sh.split('%s --help' % self.__dict__[ex]), stdout=sp.PIPE, stderr=sp.PIPE)
                p.communicate()
                self.logger.log(lg.INFO, "%s check.. OK" % ex)
            except OSError as e:
                self.logger.log(lg.CRITICAL, "could not resolve %s path: %s" % (ex, self[ex]))
                raise e

    def parse_input_file(self):
        with open(self.input_file) as SMAP:
            rdr = csv.DictReader()
            if rdr.fieldnames[0] != 'path':
                raise ValueError('first column of input csv must be a "path" column')
            self.features = [(f, set([])) for f in rdr.fieldnames[2:]]
            self.samples = {}
            for sline in rdr:
                sp = sline['path'].strip()
                if not os.path.isfile(sp):
                    raise ValueError('Could not access path %s.' % sp)
                s = Sample(**sline)
                s.order = self.features
                for f, valset in self.features: valset.add(sline[f])

    def generate_dir_tree(self):
        if self.output_dir is None:
            self.output_dir = canonic_path('./%s' % ANALYSIS_OUTPUT)
        d = self.output_dir
        self.tmp_dir = d + os.sep + TMP_NAME
        self.annotdata_dir = d + os.sep + ANNOTDATA_NAME
        self.bw_dir = d + os.sep + BIGWIG_NAME
        if self.keep_splitbam:
            self.bam_dir = d + os.sep + self.bam_dirname
        if os.path.isdir(self.tmp_dir): shutil.rmtree(self.tmp_dir)

        self.create_dir_and_log(d, lg.INFO)
        self.create_dir_and_log(self.bam_dir)
        self.create_dir_and_log(self.tmp_dir)
        self.create_dir_and_log(self.bw_dir)
        self.create_dir_and_log(self.data_dir)

    def aftermath(self):
        # remove temp folder
        # make everything read only (optional?)
        # store pipeline code?
        # merge and report statistics
        # merge data to single (usable) files

        self.logger.log(lg.CRITICAL, 'All done.')
        self.copy_log()

    def create_dir_and_log(self, path, level=lg.DEBUG):
        create_dir(path)
        self.logger.log(level, 'creating folder %s' % path)


    def parse_sample_db(self):
        feats, fmap, samples, = {}, {}, {}
        rdr = csv.DictReader(self.sample_db)
        assert rdr.fieldnames[0] == 'barcode',
        for f in rdr.fieldnames[2:]:
            m = re.match('\s+(\w+)\s+[(\w+),?(str|int|float)].*', f)
            if m is not None:
                SampleFeature(m.groups(1), units=m.groups(2), type=m.groups(3))
                feats.append()


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('experiment_path', type=str, default=None,
                   help='path to the experiment folder, is should contain a sampledb file and a BAM folder.')
    g.add_argument('--sample_db', '-sd', type=str, default=None,
                   help='path to the sample db of the experiment, defatuls to "experiment_path/sample_db.csv"')

    g = p.add_argument_group('Splitting')
    g.add_argument('--split', '-S', type=str, default='{str : {w: str()+, c: str()-, both: None}}',
                   help='A specification of in-silico splitting of each sample to sub-samples, e.g.:'
                        '{pA : {true: polyA()+, all: None},\n'
                        ' str : {w: str()+, c: str()-, both: None}}\n'
                        'will split each sample to 6 sub samples: sampl-name_pA-(y|all)_str-(w|c|both).\n'
                        'See --split_spec for available splitting options')
    g.add_argument('--split_spec', '-sh', type=str, action='store_true',
                   help='print splitting help and exit')

    g = p.add_argument_group('Reducing')
    g.add_argument('--sample_reducers', '-sr', type=str, default='',
                   help='the opposite of splitting, this allows you to ignore samples in downstream analyses, to merge'
                        'samples according to some criteria, etc. For example, the following command: "sum(X)"'
                        ' will simply sum all samples over their X feature values, and produce samples without the X '
                        'feature. sum(X=[x1,x2]) will do the same, only summing over specific values of X, while ignoring '
                        'all other values. See --reduce_spec for more info.')
    g.add_argument('--reduce_spec', '-rh', type=str, action='store_true',
                   help='print reduciton help and exit')

    g = p.add_argument_group('Project')
    g.add_argument('--project', '-P', type=str, default='cov:cov(), 3p:3p()',
                   help="This option determines how sample reads are projected into genomic functions (bp->value). Each "
                        "comma separated counter will produce a different output per sample. Projection outputs are "
                        "written as bed.gz files in the %s folder. For example 'cen11:centers(11)' will only sum the "
                        "read centers that are extended by 5 bp on each side, and will output the result to "
                        "'<sample_name>.cov11.bed.gz'. "
                        "See --project_spec for available projection options. " % PROJECTION_FOLDER)
    g.add_argument('--project_spec', '-ph', type=str, action='store_true',
                   help='print projection help and exit')

    g = p.add_argument_group('Annotation')
    g.add_argument('--annotate', '-A', type=str, default=None,
                   help="A collection of comma separated commands to generate aligned genomic arrays. For example,"
                        "{TSS(/path/to/tss.csv,[-200,200]) > mat, REB1(/path/to/tts.csv,[-300,200]) > numpy}. The given"
                        "paths are to csv files with a row per annotation, with (id, )chr, pos, strand. "
                        "For each annotation command, a single file is generated in --data_path, that contains all the "
                        "requested information in a format defined by the output type (see --annot_spec). ")
    g.add_argument('--project_spec', '-ph', type=str, action='store_true',
                   help='print projection help and exit')

    # TODO:? g = p.add_argument_group('MetaPlots')
    # g.add_argument('--meta_plot', '-M', type=str, default=None,
    #                help="A collection of comma separated commands to generate meta plots of aligned genomic features."
    #                     "For example, {tss_by_txn(/path/to/txn_tss.csv,(-200,200),sample_reducer[,clrs]]) > png}, will produce a "
    #                     "tss_by_txn meta plot. The txn_tss.csv file is as in --annotate with a grouping column - these"
    #                     "groups define the different lines plotted. The plotted surrounding genomic region is detailed "
    #                     "next, and the sample reducer needs to be defined and passed in the --sample_reduce option. This "
    #                     "reducer is applied to the samples to obtain
    #                     "where every line correspondspaths are to csv files with a row per annotation, with (id, )chr, pos, strand. "
    #                     "For each annotation command, a single file is generated in --data_path, that contains all the "
    #                     "requested information in a format defined by the output type (see --annot_spec). ")
    # g.add_argument('--project_spec', '-ph', type=str, action='store_true',
    #                help='print projection help and exit')

    g = p.add_argument_group('Browser Hubs')
    g.add_argument('--gb_hubs', '-H', type=str, default=None,
                   help="each hub specification string given will generate a hub. For example:"
                        "dko_salt(reduce=avg(repeat), normalize=q(time), collate_by=time, clr_by=dX), will generate"
                        "a hub named dko_salt, in which the samples are quantile normalized over time, the repeats are "
                        "averaged, each track shows several lines for the different time samples(from light to dark),"
                        "and the base color for each track is determined by the dX variable (i.e. first deletion)."
                        "Each hub will produce a folder of bigwig files under the '%s' folder. If a --www_path is given"
                        "a link from the wwwpath to the generated folder is added, tested, and reported."
                        "Using pre-defined reducers, " % HUB_FOLDER)
    g.add_argument('--project_spec', '-ph', type=str, action='store_true',
                   help='print projection help and exit')

    g = p.add_argument_group('Output')
    g.add_argument('--www_path', '-wp', type=str, action='store',
                   help='If given, hub linkes will be added to this path')
    g.add_argument('--user_emails', '-ue', default=None, type=str,
                   help="if provided these comma separated emails will receive notifications of ctitical "
                        "events (checkpoints, fatal errors, etc.)")
    g.add_argument('--debug', '-d', default=None, type=str,
                   help='Highly recommended. Use this mode with a pair of comma separated integers:'
                        '<numreads>,<numsamples>. The pipeline will extract this number of lines from '
                        'every input file pair, and will only run with this number of samples. If output folder'
                        '(-od) is not given, results are written to a "debug" folder.')

    g = p.add_argument_group('Execution and third party')
    g.add_argument('--exec_on', '-eo', choices=['slurm', 'local'], default='slurm',
                   help='whether to to submit tasks to a slurm cluster (default), or just use the local processors')
    g.add_argument('--samtools_exec', '-sx', type=str, default=SAMTOOLS_EXEC,
                   help='full path to samtools executable'),
    g.add_argument('--bedtools_exec', '-sx', type=str, default=SAMTOOLS_EXEC,
                   help='full path to bedtools executable'),
    g.add_argument('--bed2wig_exec', '-b2wx', type=str, default=BG2W_EXEC,
                   help='full path to badGraphToBigWig executable'),
    g.add_argument('--gencov_exec', '-gcx', type=str, default=GC_EXEC,
                   help='full path to genomeCoverageBed executable'),

    return p


def parse_args(p):
    """
    :param p: argument parser
    :return: the arguments parsed, after applying all argument logic and conversion
    """
    args = p.parse_args()
    args.__dict__['start_after'] = USER_STATES[args.start_after]
    if args.start_after != BEGIN:
        if args.output_dir is None:
            print('If the --start_after option is used, an existing output directory must be provided (-od).')
            exit()
            args.__dict__['fastq_pref'] = args.output_dir  # ignoring input folder
    else:
        if args.fastq_prefix is None:
            print(
                'If the --start_after option is not used, an input fastq prefix/folder mut be provided (--fastq_prefix).')
            exit()

    p, s = os.path.split(args.fastq_prefix)
    args.__dict__['fastq_path'] = p
    args.__dict__['fastq_pref'] = s
    args.__dict__['fastq_prefix'] = canonic_path(args.fastq_prefix)
    p, s = os.path.split(args.fastq_prefix)
    if args.barcode_file is None:
        args.__dict__['barcode_file'] = os.sep.join([args.fastq_path, 'barcodes'])

    if args.debug is not None:
        nlines, nsamples = args.debug.split(',')
        args.__dict__['db_nlines'] = int(nlines)
        args.__dict__['db_nsamples'] = int(nsamples)

    if args.output_dir is not None:
        args.__dict__['output_dir'] = canonic_path(args.__dict__['output_dir'])

    args.__dict__['filters'] = build_filter_schemes(args.filters)

    analyzers, sel_a = collect_analyzers(), {}
    for a in args.analyses.split(','):
        if a not in analyzers:
            print("could not resolve analyzer '%s', use run -ah for more details" % a)
            exit()
        else: sel_a[a] = analyzers[a]  # still need to instantiate them
    args.__dict__['analyzers'] = sel_a

    return args


if __name__ == '__main__':
    p = build_parser()
    a = parse_args(p)
    mh = MainHandler(a, ' '.join(sys.argv))
    mh.execute()