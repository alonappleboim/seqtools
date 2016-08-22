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
from singleend_filters import *


BEGIN = 0
FASTQ = 1
ALIGN = 2
FILTR = 3
SPLIT = 4

STATES = ['BEGIN', 'FASTQ', 'ALIGN', 'FILTR', 'SPLIT']
USER_STATES = {'BEGIN':BEGIN, 'FASTQ':FASTQ, 'ALIGN':ALIGN, 'FILTR':FILTR, 'SPLIT':SPLIT}


class MainHandler(object):

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

    def __init__(self, argobj, cmdline):
        self.__dict__.update(argobj.__dict__)

        self.setup_log()
        self.comq = mp.Queue()
        if self.debug:
            self.logger.log(lg.INFO, '=================== DEBUG MODE (%s) ===================' % self.debug)
        self.check_third_party()

        msg = 'filters:\n' + '\n'.join(scheme.name+'\n'+str(scheme) for scheme in self.filters.values())
        self.logger.log(lg.INFO, msg)

        self.parse_barcode_file()
        self.generate_dir_tree()
        shutil.copy(self.barcode_file, self.output_dir + os.sep + 'barcodes')
        self.files = {S: {s: {} for s in self.s2b.keys()}
                      for S in USER_STATES.values()}

        instantiated = {}
        for a, acls in self.analyzers.items():
            instantiated[a] = acls(self.tmp_dir)
        self.analyzers = instantiated
        msg = 'analyzers:\n' + '\n'.join(a for a in self.analyzers.keys())
        self.logger.log(lg.INFO, msg)

        self.w_manager = WorkManager(self.comq, self.exec_on=='slurm')

    def execute(self):
        if self.start_after <= BEGIN: self.make_fastq()
        if self.start_after <= FASTQ: self.make_bam()
        if self.start_after <= ALIGN: self.filter()
        if self.start_after <= FILTR: self.analyze()
        self.aftermath()

    def make_fastq(self):
        self.logger.log(lg.INFO, 'Re-compiling fastq files...')
        fs = self.files[FASTQ]
        self.collect_input_fastqs()
        bcout = None
        if self.keep_nobarcode:
            bcout = self.fastq_dir + os.sep + NO_BC_NAME + '.fastq.gz'
        self.split_barcodes(no_bc=bcout)
        token_map = {}
        for bc, sample in self.b2s.items():

            fs[sample]['in1'] = self.tmp_dir + os.sep + sample + '-1'
            fs[sample]['in2'] = self.tmp_dir + os.sep + sample + '-2'
            fs[sample]['fastq'] = self.fastq_dir + os.sep + sample + '.fastq.gz'
            args = dict(files=fs[sample], bc_len=self.bc_len, umi_len=self.umi_length)
            token_map[self.w_manager.run(format_fastq, kwargs=args)] = sample
        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error in processing sample %s:\n%s' % (sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                # in principal the next step could start now, but then recovery options make the code
                # practically unreadable, and the performance benefit is very small
                self.logger.log(lg.DEBUG, '%s ready.' % fs[sample]['fastq'])
        self.logger.log(lg.INFO, 'fastq files were written to: %s' % self.fastq_dir)
        self.mark(FASTQ)

    def make_bam(self):
        self.logger.log(lg.INFO, 'Aligning to genome...')
        token_map = {}
        for bc, sample in self.b2s.items():
            f = self.files[ALIGN][sample]
            f['fastq'] = self.files[FASTQ][sample]['fastq']
            f['bam'] = self.bam_dir + os.sep + sample + '.bam'
            f['tmp_bam'] = self.tmp_dir + os.sep + sample + TMP_BAM_SUFF
            f['sam_hdr'] = self.bam_dir + os.sep + sample + SAM_HDR_SUFF
            f['align_stats'] = self.tmp_dir + os.sep + sample + BT_STATS_SUFF
            if self.klac_index_path is not None:
                f['klac_align_stats'] = self.tmp_dir + os.sep + sample + '.klac' + BT_STATS_SUFF
            if self.keep_unaligned:
                f['unaligned_bam'] = self.unaligned_dir + os.sep + sample + '.bam'
            args = dict(files=self.files[ALIGN][sample], bowtie_exec=self.bowtie_exec,
                        n_threads=self.n_threads, scer=self.scer_index_path, klac=self.klac_index_path)
            token_map[self.w_manager.run(align, kwargs=args)] = sample

        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error in aligning sample %s:\n%s' % (sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                # in principal the next step could start now, but then recovery options make the code
                # practically unreadable, and the performance benefit is very small
                self.logger.log(lg.DEBUG, '%s ready.' % self.files[ALIGN][sample]['bam'])

        for f in os.listdir(self.tmp_dir):
            if not f.endswith(BT_STATS_SUFF): continue
            is_klac = f.endswith('.klac'+BT_STATS_SUFF)
            sample = f.split('.')[0]
            fpath = self.tmp_dir + os.sep + f
            with open(fpath) as F:
                for line in F:
                    stat, cnt = line.strip().split('\t')
                    stat += '_klac' if is_klac else ''
                    if stat not in self.stat_order: self.stat_order.append(stat)
                    self.stats[sample][stat] += int(cnt)
            os.remove(fpath)
        self.mark(ALIGN)

    def filter(self):
        self.logger.log(lg.INFO, 'Filtering the reads...')
        token_map = {}
        for bc, sample in self.b2s.items():
            for f, fscheme in self.filters.items():
                files = self.files[FILTR][sample][f] = {}
                files['bam_in'] = self.files[ALIGN][sample]['bam']
                files['sam_hdr'] = self.files[ALIGN][sample]['sam_hdr']
                files['bam_out'] = os.sep.join([self.output_dir, f, self.bam_dirname, sample + '.bam'])
                args = dict(files=files, fscheme=fscheme)
                token_map[self.w_manager.run(filter_bam, kwargs=args)] = (sample, f)

        while token_map:
            tok, e, info = self.comq.get()
            sample, filter = token_map.pop(tok)
            if e is not None:
                msg = 'error while filtering sample %s with filter %s:\n%s' % (sample, filter, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                # in principal the next step could start now, but then recovery options make the code
                # practically unreadable, and the performance benefit is very small
                stat = 'passed:%s' % f
                if stat not in self.stat_order: self.stat_order.append(stat)
                self.stats[sample][stat] += info
                self.logger.log(lg.DEBUG, 'filtered %s with %s (%i passed).' % (sample, f, info))

        self.mark(FILTR)

    def analyze(self):
        self.logger.log(lg.INFO, 'Analyzing the data...')
        token_map = {}
        for bc, sample in self.b2s.items():
            for f, fscheme in self.filters.items():
                for aname, analyzer in self.analyzers.items():
                    files = self.files[FILTR][sample][f][aname] = {}
                    files['bam_in'] = self.files[FILTR][sample][f]['bam_out']
                    for s in analyzer.suffixes.keys():
                        files[s] = os.sep.join([self.output_dir, analyzer.out_folder, sample + '_' + f + s])
                    args = dict(files=files, analyzer=analyzer)
                    token_map[self.w_manager.run(analyze_bam, kwargs=args)] = (sample, f, aname)
        while token_map:
            tok, e, info = self.comq.get()
            sample, filter, analyzer = token_map.pop(tok)
            if e is not None:
                msg = 'error while analyzing (%s, %s) with %s analysis:\n %s' % (sample, filter, analyzer, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                # in principal the next step could start now, but then recovery options make the code
                # practically unreadable, and the performance benefit is very small
                if info:
                    for s, cnt in info.items():
                        stat = '%s:%s' % (filter, s)
                        if stat not in self.stat_order: self.stat_order.append(stat)
                        self.stats[sample][stat] += info
                self.logger.log(lg.DEBUG, 'analyzed (%s, %s) with %s.' % (sample, f, analyzer))

    def print_stats(self):
        stats = set([])
        fns = ['sample'] + self.stat_order
        fh = open(self.output_dir + os.sep + self.exp + '_statistics.csv','w')
        wrtr = csv.DictWriter(fh, fns)
        wrtr.writeheader()
        for s, st in self.stats.items():
            wrtr.writerow(dict(sample=s, **{k:str(st[k]) for k in fns[1:]}))
        fh.close()

    def mark(self, stage):
        self.print_stats()
        stage = STATES[stage]
        msg = ('Finished stage %s. You can continue the pipeline from this point '
               'with the option -sa %s (--start_after %s)' % (stage, stage, stage))
        self.logger.log(lg.INFO, msg)
        self.copy_log()
        fh = open(self.output_dir + os.sep + '.pipeline_state', 'w')
        fh.write('%s\n' % stage)
        fh.close()

    def get_mark(self):
        fh = open(self.output_dir + os.sep + '.pipeline_state')
        stage = int(fh.read().strip())
        fh.close()
        return stage

    def copy_log(self):
        shutil.copy(self.logfile, self.output_dir + os.sep + 'full.log')

    def prepare_outputs(self):
        pass

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

    def parse_barcode_file(self):
        with open(self.barcode_file) as BCF:
            exp = re.match('.*experiment.*:\s+(\w+)', BCF.readline())
            if exp is None:
                msg = 'barcodes file should contain a header with experiment name: ' \
                      'experiment: <expname>'
                self.logger.log(lg.CRITICAL, msg)
                raise ValueError(msg)
            self.user = getpass.getuser()
            self.exp = exp.group(1)
            msg = 'user: %s, experiment: %s' % (self.user, self.exp)
            self.logger.log(lg.INFO, msg)
            b2s, s2b = {}, {}
            bl = None
            for i, line in enumerate(BCF):
                if self.debug is not None:
                    if i >= self.db_nsamples: break #limit number of samples
                b, s = line.strip().split(',')
                if b in b2s:
                    msg = 'barcode %s is not unique.' % b
                    self.log(msg, lg.CRITICAL)
                    raise(ValueError(msg))
                if s in s2b:
                    msg = 'sample %s is not unique.' % s
                    self.log(msg, lg.CRITICAL)
                    raise (ValueError(msg))
                if bl is None: bl = len(b)
                elif bl != len(b):
                    msg = 'barcode %s does not conform to barcode length %i.' % (b, bl)
                    self.log(msg, lg.CRITICAL)
                    raise (ValueError(msg))
                b2s[b] = s
                s2b[s] = b

            msg = '\n'.join(['barcodes:'] + ['%s -> %s' % bs for bs in b2s.items()])
            self.logger.log(lg.DEBUG, msg)
            self.logger.log(lg.INFO, 'found %i samples.' % len(b2s))

        self.b2s = b2s
        self.s2b = s2b
        self.bc_len = len(list(b2s)[0])
        self.stats = {s: Counter() for s in s2b.keys()}
        self.stat_order = []

    def write_barcode_file(self):
        shutil.copy(self.barcode_file, self.output_dir + os.sep + 'barcodes')

    def collect_input_fastqs(self):
        files = {}

        for fn in os.listdir(self.fastq_path):
            if not fn.endswith('fastq.gz'): continue
            if not fn.startswith(self.fastq_pref): continue
            parts = re.split('_R\d',fn)
            if len(parts) == 2:
                path = self.fastq_path + os.sep + fn
                pref = parts[0]
                if pref in files:
                    if '_R1' in files[pref]:
                        files[pref] = (files[pref], path)
                    else:
                        files[pref] = (path, files[pref])
                else:
                    files[pref] = path

        files = [f for f in files.values() if type(()) == type(f)]
        msg = '\n'.join('found fastq files:\n%s\n%s' % fs for fs in files)
        self.logger.log(lg.INFO, msg)
        if not files:
            msg = "could not find R1/R2 fastq.gz pairs in given folder: %s" % self.fastq_path
            self.logger.log(lg.CRITICAL, msg)
            raise IOError(msg)
        self.input_files = files

    def generate_dir_tree(self):
        if self.start_after != BEGIN:
            try:
                cur = self.get_mark()
                if USER_STATES[cur] >= USER_STATES[self.start_after]:
                    msg = 'restarting from %s in folder: %s ' % (self.start_after, self.output_dir)
                    self.logger.log(lg.INFO, msg)
                else:
                    msg = 'folder state %s in folder %s incompatible with start_after %s request' \
                          % (cur, self.output_dir, self.start_after)
                    self.logger.log(lg.CRITICAL, msg)
                    exit()
            except IOError:
                msg = 'could not find an existing output folder: %s' % self.output_dir
                self.logger.log(lg.CRITICAL, msg)
                exit()
        else:
            if self.output_dir is None:
                folder = canonic_path(DATA_PATH) + os.sep + self.user
                create_dir(folder)
                folder += os.sep + self.exp
                create_dir(folder)
                folder += os.sep + datetime.datetime.now().strftime("%d-%m-%y")
                create_dir(folder)
                if self.debug is not None:
                    fname = 'debug'
                else:
                    i = [int(x) for x in os.listdir(folder) if isint(x)]
                    if not i: i = 1
                    else: i = max(i) + 1
                    fname = str(i)
                folder += os.sep + fname
                self.output_dir = folder

        d = self.output_dir
        self.tmp_dir = d + os.sep + TMP_NAME
        self.fastq_dir = d + os.sep + self.fastq_dirname
        self.bam_dir = d + os.sep + self.bam_dirname
        if os.path.isdir(self.tmp_dir): shutil.rmtree(self.tmp_dir)

        if self.start_after == BEGIN:
            # assuming all folder structure exists if check passes
            self.create_dir_and_log(d, lg.INFO)
            self.create_dir_and_log(self.fastq_dir)
            self.create_dir_and_log(self.bam_dir)
            self.create_dir_and_log(self.tmp_dir)
            if self.keep_unaligned:
                self.unaligned_dir = d + os.sep + UNALIGNED_NAME
                self.create_dir_and_log(self.unaligned_dir)

        if self.start_after <= ALIGN:
            for f, scheme in self.filters.items():
                create_dir(os.sep.join([self.output_dir, f]))
                filter_folder = os.sep.join([self.output_dir, f, self.bam_dirname])
                create_dir(filter_folder)
                scheme.folder = filter_folder

        for a, acls in self.analyzers.items():
            create_dir(os.sep.join([self.output_dir, acls.out_folder]))

    def split_barcodes(self, no_bc=None):
        #
        #  first - compile awk scripts that will split input according to barcodes,
        # then use them with some smart pasting. For each input pair:
        # paste <(zcat {R1}) <(zcat {R2}) |\ %paste both files one next to the other
        # paste - - - -|\ % make every 4 lines into one
        # awk -F "\\t" -f {script1-exact_match}' |\ % split exact barcode matches
        # awk -F "\\t" -f {script2-ham_dist}' % split erroneous matches (subject to given hamming distance)
        #
        """

        :param no_bc: if given, orphan fastq enries are written to this prefix (with R1/R2 interleaved)
        """
        def compile_awk(it, b2s):
            cnt_path = self.tmp_dir + os.sep + (BC_COUNTS_FNAME %  it)
            nobc = NO_BC_NAME + '-' + it
            arraydef = ';\n'.join('a["%s"]="%s-%s"' % (b, s, it) for b, s in b2s.items()) + ';\n'
            awk_str = (""" 'BEGIN {%s} {x=substr($4,1,%i); if (x in a) """,
                       """{c[a[x]]++; print >> "%s/"a[x];} else {c["%s"]++; print;} }""",
                       """END { for (bc in c) print bc, c[bc] >> "%s" } '""")
            awk_str = ''.join(awk_str) % (arraydef, self.bc_len , self.tmp_dir, nobc, cnt_path)
            return awk_str, cnt_path

        def merge_statistics(bc1, bc2):
            stat = "#reads"
            self.stat_order.append(stat)
            self.stats[NO_BC_NAME] = Counter()
            with open(bc1) as IN:
                for line in IN:
                    sample, cnt = line.strip().split(' ')
                    sample = sample[:-2]
                    if sample == 'no-barcode': continue  # only from bc_counts-2
                    self.stats[sample][stat] += int(cnt)
            os.remove(bc1)
            with open(bc2) as IN:
                for line in IN:
                    sample, cnt = line.strip().split(' ')
                    self.stats[sample[:-2]][stat] += int(cnt)
            os.remove(bc2)
            for s in self.s2b.keys():
                if s not in self.stats:
                    self.stats[s][stat] += 0
            msg = '\n'.join(['%s: %i' % (s, c['#reads']) for s, c in self.stats.items()])
            self.logger.log(lg.CRITICAL, 'read counts:\n' + msg)

        hb = {}
        for b,s in self.b2s.items():
            hb.update({eb:s for eb in hamming_ball(b, self.hamming_distance)})

        awk1p, cnt1 = compile_awk("1", self.b2s)
        awk2p, cnt2 = compile_awk("2", hb)
        outf = open(os.devnull, 'w') if no_bc is None else open(no_bc, 'wb')
        for r1, r2 in self.input_files:
            msg = 'splitting files:\n%s\n%s' % (os.path.split(r1)[1],os.path.split(r2)[1])
            self.logger.log(lg.INFO, msg)
            paste1 = sp.Popen('paste <(zcat %s) <(zcat %s)' % (r1,r2), stdout=sp.PIPE,
                              shell=True, executable='/bin/bash')
            awkin = sp.Popen(sh.split('paste - - - -'), stdin=paste1.stdout, stdout=sp.PIPE)
            if self.debug: # only a small subset of reads
                awkin = sp.Popen(sh.split('head -%i' % self.db_nlines), stdin=awkin.stdout, stdout=sp.PIPE)
            awk1 = sp.Popen(sh.split('awk -F "\\t" ' + awk1p), stdin=awkin.stdout, stdout=sp.PIPE)
            awk2 = sp.Popen(sh.split('awk -F "\\t" ' + awk2p), stdin=awk1.stdout, stdout=sp.PIPE)
            awkcmd = """awk -F "\\t" '{print $1"\\n"$3"\\n"$5"\\n"$7; print $2"\\n"4"\\n"$6"\\n"$8;}' """
            wfastq = sp.Popen(sh.split(awkcmd), stdin=awk2.stdout, stdout=sp.PIPE)
            gzip = sp.Popen(['gzip'], stdin=wfastq.stdout, stdout=outf)
        gzip.wait()
        self.logger.log(lg.INFO, 'Barcode splitting finished.')

        merge_statistics(cnt1, cnt2)

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


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('--fastq_prefix', '-fp', type=str, default=None,
                   help='path to a prefix of fastq files (R1 & R2) containing the transeq data.'
                        'This can be a folder (must end with "/"), in which case all R1/R2 pairs'
                        'in the folder are considered, or a "path/to/files/prefix", in which case '
                        'all files in the path with the prefix are considered')
    g.add_argument('--start_after', '-sa', default='BEGIN',
                   choices=[k for k in USER_STATES.keys()],
                   help='If given the pipeline will try to continue a previous run, specified through '
                        'the "output_dir" argument, from the selected stage. In this case --fastq_prefix'
                        ' is ignored.')

    g = p.add_argument_group('Output')
    g.add_argument('--output_dir', '-od', default=None, type=str,
                   help='path to the folder in which most files are written'
                        'If not given, the date and info from barcode file are used to '
                        'generate a new folder in %s' % DATA_PATH)
    g.add_argument('--fastq_dirname', '-fd', default='FASTQ', type=str,
                   help='name of folder in which fastq files are written. relative to output_dir')
    g.add_argument('--bam_dirname', '-bd', default='BAM', type=str,
                   help='name of folder in which bam files are written relative to output_dir'
                        '(or each filter folder)')
    g.add_argument('--user_emails', '-ue', default=None, type=str,
                   help="if provided these comma separated emails will receive notifications of ctitical "
                        "events (checkpoints, fatal errors, etc.)")
    g.add_argument('--debug', '-d', default=None, type=str,
                   help='Highly recommended. Use this mode with a par of comma separated integers:'
                        '<numlines>,<numsamples>. The pipeline will extract this number of lines from '
                        'every input file pair, and will only run with this number of samples out of the '
                        'given barcode file. If output folder (-od) is not given, results are written to '
                        'a "debug" folder.')

    g = p.add_argument_group('Barcode Splitting')
    g.add_argument('--barcode_file', '-bf', type=str, default=None,
                   help='a file with 2 header lines:'
                        '\nname: <name> \n experiment: <expname>\n followed by any number '
                        'of lines specifying the barcodes and samples: <bc1>,<sample1>\n<bc2>,<sample2>...\n'
                        'default is "barcodes" in given fastq path')
    g.add_argument('--umi_length', '-ul', type=int, default=8,
                   help='UMI length')
    g.add_argument('--hamming_distance', '-hd', default=1, type=int,
                   help='barcode upto this hamming distance from given barcodes are handled by '
                        'the pipeline')
    g.add_argument('--keep_nobarcode', '-kbc', action='store_true',
                   help='keep fastq entries that did not match any barcode')

    g = p.add_argument_group('Alignment')
    g.add_argument('--scer_index_path', '--sip', type=str, default='/cs/wetlab/genomics/bowtie_index/scer/sacCer3',
                   help='path prefix of s. cervisae genome bowtie index')
    g.add_argument('--klac_index_path', '-kip', type=str, default=None,
                   help='If given, data is also aligned to k. lacis genome (can be found in '
                        '/cs/wetlab/genomics/bowtie_index/klac/)')
    g.add_argument('--n_threads', '-an', type=int, default=4,
                   help='number of threads used for alignment per bowtie instance')
    g.add_argument('--keep_unaligned', '-ku', action='store_true',
                   help='if set, unaligned reads are written to '
                        'output_folder/%s/<sample_name>.bam' % UNALIGNED_NAME)

    g = p.add_argument_group('Filters and Splitters',
                             description='different filters applied to base BAM file. Each filter result is '
                                         'processeed downstream and reported separately')
    g.add_argument('--filters', '-F', action='store',
                   default='unique:dup(),qual(qmin=5)+',
                   help='specify filter schemes to apply to data. Expected string conforms to ([] are for grouping):\n' \
                        '[<filter_scheme>:<filter>([<argname1=argval1>,]+)[+|-]);]*\n the filter_scheme will be used to name all '
                        'resulting outputs from this branch of the data. use "run -fh" for more info.')
    g.add_argument('--filter_specs', '-fh', action='store_true',
                   default='print available filter specs and filter help and exit')

    g = p.add_argument_group('Filters',
                             description='different filters applied to base BAM file. Each filter result is '
                                         'processeed downstream and reported separately')
    g.add_argument('--split', '-S', action='store',
                   default='[pA, [+:polyA()+, -:polyA()-)]],[str, [w:strand()+,c:strand()-]]',
                   help='specify filter schemes to apply to data. Expected string conforms to ([] are for grouping):\n' \
                        '[<filter_scheme>:<filter>([<argname1=argval1>,]+)[+|-]);]*\n the filter_scheme will be used to name all '
                        'resulting outputs from this branch of the data. use "run -fh" for more info.')

    g = p.add_argument_group('Analyses',
                             description='different analyses applied to filtered BAM files. Each analysis is '
                                         'performed on single samples, and when all is done, results are merged '
                                         'to a single file in the output/path/results folder')
    g.add_argument('--analyses', '-A', action='store', default='3pT,covT',
                   help='specify analyses to apply to data. a comma separated list (see run -ah" for more info).')
    g.add_argument('--analysis_specs', '-ah', action='store_true',
               default='print available analyses specs and analysis help and exit')

    g = p.add_argument_group('Execution and third party')
    g.add_argument('--exec_on', '-eo', choices=['slurm', 'local'], default='slurm',
                   help='whether to to submit tasks to a slurm cluster (default), or just use the local processors')
    g.add_argument('--bowtie_exec', '-bx', type=str, default=BOWTIE_EXEC,
                   help='full path to bowtie executable'),
    g.add_argument('--samtools_exec', '-sx', type=str, default=SAMTOOLS_EXEC,
                   help='full path to samtools executable'),
    g.add_argument('--bed2wig_exec', '-b2wx', type=str, default=BG2W_EXEC,
                   help='full path to badGraphToBigWig executable'),
    g.add_argument('--gencov_exec', '-gcx', type=str, default=GC_EXEC,
                   help='full path to genomeCoverageBed executable'),

    g = p.add_argument_group('polyA')
    g.add_argument('--polyA', '-pa', action='store_false',
                   help='set to avoid polyA reads special treatment'),
    g.add_argument('--polyA_threshold', '-pt', type=int, default=4,
                   help='any read with more than this number of As at its *unaligned* end is considered a polyA read')
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