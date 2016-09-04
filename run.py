
import stat
from collections import Counter, OrderedDict
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
import pickle

from config import *

from workers import *
from exporters import *
from utils import *
from filters import *

if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    import subprocess as sp
    import os
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()


BEGIN = 0
FASTQ = 1
BAM = 2

STATES = ['BEGIN', 'FASTQ', 'BAM']
USER_STATES = {'BEGIN':BEGIN, 'FASTQ':FASTQ, 'BAM':BAM}


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
        self.args = argobj

        self.setup_log()
        self.logger.log(lg.INFO, 'commandline: %s' % cmdline)
        self.comq = mp.Queue()
        if self.debug:
            self.logger.log(lg.INFO, '=================== DEBUG MODE (%s) ===================' % self.debug)
        self.check_third_party()

        self.bc_len, self.samples, self.features = self.parse_sample_db()
        self.generate_dir_tree()
        sfname = self.output_dir + os.sep + 'sample_db.csv'
        if not os.path.isfile(sfname): shutil.copy(self.sample_db, sfname)

        self.stats = OrderedDict()
        for s in self.samples.values(): self.stats[s.base_name()] = Counter()
        self.stat_order = []
        self.fpipe = build_filter_schemes('filter:'+self.filter)['filter']
        self.logger.log(lg.INFO, 'Filters:\n' + str(self.fpipe))
        self.exporters = exporters_from_string(self.exporters, self.output_dir)
        self.w_manager = WorkManager(self.comq, self.exec_on=='slurm')

    def execute(self):
        if self.start_after <= BEGIN: self.make_fastq()
        if self.start_after <= FASTQ: self.make_bam()
        self.track()
        self.count()
        self.export()
        self.aftermath()

    def make_fastq(self):
        self.logger.log(lg.INFO, 'Re-compiling fastq files...')
        self.collect_input_fastqs()
        bcout = None
        if self.keep_nobarcode:
            bcout = self.fastq_dir + os.sep + NO_BC_NAME + '.fastq.gz'
        self.split_barcodes(no_bc=bcout)
        token_map = {}
        for bc, sample in self.samples.items():
            sf = sample.files
            sf['in1'] = self.tmp_dir + os.sep + sample.base_name() + '-1'
            sf['in2'] = self.tmp_dir + os.sep + sample.base_name()+ '-2'
            sf['fastq'] = self.fastq_dir + os.sep + sample.base_name() + '.fastq.gz'
            args = dict(files=sf, bc_len=self.bc_len, umi_len=self.umi_length)
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
                self.logger.log(lg.DEBUG, '%s ready.' % sf['fastq'])
        self.logger.log(lg.INFO, 'fastq files were written to: %s' % self.fastq_dir)
        self.checkpoint(FASTQ)

    def make_bam(self):
        self.logger.log(lg.INFO, 'Aligning to genome...')
        token_map = {}
        for bc, sample in self.samples.items():
            sf = sample.files
            sf['bam'] = self.bam_dir + os.sep + sample.base_name() + '.bam'
            sf['unfiltered_bam'] = self.bam_dir + os.sep + sample.base_name() + '.unfiltered.bam'
            sf['bam_f'] = None
            if self.keep_filtered:
                sf['bam_f'] = self.filtered_dir + os.sep + sample.base_name() + '.bam'
            sf['tmp_bam'] = self.tmp_dir + os.sep + sample.base_name() + TMP_BAM_SUFF
            sf['sam_hdr'] = self.bam_dir + os.sep + sample.base_name() + SAM_HDR_SUFF
            sf['align_stats'] = self.tmp_dir + os.sep + sample.base_name() + BT_STATS_SUFF
            if self.klac_index_path is not None:
                sf['klac_align_stats'] = self.tmp_dir + os.sep + sample.base_name() + '.klac' + BT_STATS_SUFF
            if self.keep_unaligned:
                sf['unaligned_bam'] = self.unaligned_dir + os.sep + sample.base_name() + '.bam'
            args = dict(files=sf, bowtie_exec=self.bowtie_exec, fpipe=self.fpipe,
                        n_threads=self.n_threads, scer=self.scer_index_path, klac=self.klac_index_path)
            token_map[self.w_manager.run(make_bam, kwargs=args)] = sample

        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error in aligning sample %s:\n%s' % (sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                # in principal the next step could start now, but then recovery options make the code
                # practically unreadable, and the performance benefit is very small
                if 'passed_filter' not in self.stat_order:
                    self.stat_order.append('passed_filter')
                self.stats[sample.base_name()]['passed_filter'] = info
                self.logger.log(lg.DEBUG, '%s ready.' % sf['bam'])

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
        self.checkpoint(BAM)

    def track(self):
        # TODO: this will be replaced by a proper analysis tool, for now it's just here to support current functionality
        self.logger.log(lg.INFO, 'Making tracks...')
        token_map = {}
        for bc, sample in self.samples.items():
            sf = sample.files
            sf['cbw'] = self.bw_dir + os.sep + sample.base_name() + '.c.bw'
            sf['wbw'] = self.bw_dir + os.sep + sample.base_name() + '.w.bw'
            sf['tmp_bed'] = self.tmp_dir + os.sep + sample.base_name() + '.tmp.bed'
            token_map[self.w_manager.run(make_tracks, kwargs=dict(files=sf))] = sample

        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error in making tracks for sample %s:\n%s' % (sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                self.logger.log(lg.DEBUG, '%s ready.' % sf['cbw'])
                self.logger.log(lg.DEBUG, '%s ready.' % sf['wbw'])

        if self.make_hub: self.build_hub()

    def count(self):
        pass
        # TODO: this will be replaced by a proper analysis tool, for now it's just here to support current functionality
        self.logger.log(lg.INFO, 'Counting reads...')

        # prepare input bed file for all workers
        tmp_bed = self.tmp_dir + os.sep + 'tts.tmp.bed'
        tts_bed = open(tmp_bed, 'w')
        chrlens = chr_lengths()
        self.tts_accs = []
        for line in open(self.tts_file):
            acc, chr, orf_start, orf_end, tts = line.strip().split('\t')
            if chr == 'chrXVII': chr = 'chrM'
            orf_start, orf_end = int(orf_start), int(orf_end)
            strand = orf_start < orf_end
            tts_annot_problem = tts == 'NaN' or \
                                (strand and int(tts) < orf_end) or \
                                (not strand and int(tts) > orf_end)
            tts = orf_end if tts_annot_problem else int(tts)
            w_fr = max(0, tts + ((-1) ** (1 - strand)) * self.count_window[1 - strand])
            w_to = min(chrlens[chr], tts + ((-1) ** (1-strand)) * self.count_window[strand])
            if not self.dont_bound_start:
                if strand: w_fr = max(w_fr, orf_start)
                else: w_to = min(w_to, orf_start)
            tts_bed.write('\t'.join([chr, str(w_fr), str(w_to), acc, '1', '+' if strand else '-'])+'\n')
            self.tts_accs.append(acc)
        tts_bed.close()

        token_map = {}
        for bc, sample in self.samples.items():
            sf = sample.files
            sf['tmp_cnt'] = self.tmp_dir + os.sep + sample.base_name() + '.tmp.cnt'
            args = dict(annot_file=tmp_bed, files=sf)
            token_map[self.w_manager.run(count, kwargs=args)] = sample

        cnts = {}
        while token_map:
            tok, e, info = self.comq.get()
            sample = token_map.pop(tok)
            if e is not None:
                msg = 'error in counting sample %s:\n%s' % (sample, e)
                self.logger.log(lg.CRITICAL, msg)
            else:
                cnts[sample] = info
                self.logger.log(lg.DEBUG, 'Collected counts for sample %s' % sample)

        if self.keep_tts_bed:
            shutil.copy(tmp_bed, self.output_dir)
        else:
            os.remove(tmp_bed)

        self.tts_cnts = OrderedDict()
        for s in self.samples.values(): self.tts_cnts[s] = cnts[s] # maintaining sample_db order

    def export(self):
        stats = OrderedDict()
        for s in self.samples.values():
            stats[s] = self.stats[s.base_name()]
        s = Sample()
        for f in self.features.values(): s.fvals[f] = NO_BC_NAME
        stats[s] = self.stats[NO_BC_NAME]
        all = [('stats', stats, self.stat_order),
               ('tts', self.tts_cnts, self.tts_accs)]
        for e in self.exporters:
            fs = e.export(self.features.values(), self.samples.values(), all)
            for f in fs:
                self.logger.log(lg.INFO, 'Exported data to file: %s' % (self.output_dir+os.sep+f,))
                if self.export_path is not None:
                    target = self.export_path + os.sep + self.exp + '-' + f
                    shutil.copy(self.output_dir + os.sep + f, target)
                    self.logger.log(lg.DEBUG, 'Copied data to: %s' % target)

    def print_stats(self):
        stats = set([])
        fns = ['sample'] + self.stat_order
        fh = open(self.output_dir + os.sep + self.exp + '_statistics.csv','w')
        wrtr = csv.DictWriter(fh, fns)
        wrtr.writeheader()
        for s, st in self.stats.items():
            wrtr.writerow(dict(sample=s, **{k:str(st[k]) for k in fns[1:]}))
        fh.close()

    def checkpoint(self, stage):
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
        with open(self.output_dir + os.sep + '.pipeline_state') as fh:
            return fh.read().strip()

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

    def parse_sample_db(self):

        def parse_features(hdr):
            feat_pat = re.compile('\s*(?P<name>\w+)\s*(?:\((?P<short_name>\w+)\))?'
                                  '\s*:(?P<type>\w+)(:?\[(?P<units>\w+)\])?')
            hdr = hdr.split(DELIM)
            features = FeatureCollection()
            f_pos_map = {}
            for i, f in enumerate(hdr):
                if i == 0:
                    assert hdr[0] == 'barcode', 'first column in sample db needs to be the "barcode" column'
                elif f.startswith('#'):
                    msg = ("ignoring column %s in sample db" % f)
                    self.logger.log(lg.INFO, msg)
                else:
                    m = re.match(feat_pat, f)
                    if m is None:
                        msg = ("couldn't understand feature '%s' in sample_db file, format should be: "
                               "<name>(<short_name>):(str|int|float)[units] (short_name and units are optional) or "
                               "column is ignored if it starts with '#'" % f)
                        self.logger.log(lg.CRITICAL, msg)
                        raise (ValueError(msg))
                    try:
                        f_pos_map[i] = Feature(**m.groupdict())
                    except ValueError:
                        snames = '\n'.join(f.short_name for f in features.values)
                        msg = ("features must have distinct names and short_names - %s appears at least twice (or "
                               "its short_name matched a previous generated short name):\n%s" % f, snames)
                        self.logger.log(lg.CRITICAL, msg)
                        raise (ValueError(msg))
                    features.add_feature(f_pos_map[i])
            return features, f_pos_map

        def parse_samples(file, f_pos_map):
            b2s, bc_len = OrderedDict(), None
            for i, line in enumerate(file):
                if line.strip()[0] == '#': continue  # comment
                if self.debug is not None:
                    if i >= self.db_nsamples: break  # limit number of samples
                sample = Sample()
                for j, val in enumerate(line.strip().split(DELIM)):
                    val = val.strip()
                    if j == 0:
                        if i == 0:
                            bc_len = len(val)  # first barcode
                        elif bc_len != len(val):
                            msg = "barcode %s has a different length" % val
                            self.logger.log(lg.CRITICAL, msg)
                            raise (TypeError(msg))
                        if val in b2s:
                            msg = "barcode %s is not unique" % val
                            self.logger.log(lg.CRITICAL, msg)
                            raise (TypeError(msg))
                        sample.barcode = val
                    elif j in f_pos_map:
                        f = f_pos_map[j]
                        try:
                            v = f.type(val)
                        except ValueError:
                            msg = ("couldn't cast value %s in sample %i, feature '%s' to "
                                   "given type - %s." % (val, i + 1, f.name, f.strtype))
                            self.logger.log(lg.CRITICAL, msg)
                            raise (ValueError(msg))
                        f.vals.add(v)
                        sample.fvals[f] = v
                    if hash(sample) in [hash(s) for s in b2s.values()]:
                        msg = "2 samples (or more) seem to be identical - %s" % sample
                        self.logger.log(lg.CRITICAL, msg)
                        raise (TypeError(msg))
                b2s[sample.barcode] = sample
            return b2s, bc_len

        sdb = open(self.sample_db)
        exp = re.match('.*experiment.*:\s+(\w+)', sdb.readline())
        if exp is None:
            msg = 'barcodes file should contain a header with experiment name: ' \
                  'experiment: <expname>'
            self.logger.log(lg.CRITICAL, msg)
            raise ValueError(msg)
        self.user = getpass.getuser()
        self.exp = exp.group(1)
        msg = 'user: %s, experiment: %s' % (self.user, self.exp)
        self.logger.log(lg.INFO, msg)
        features, f_pos_map = parse_features(sdb.readline())
        b2s, bc_len = parse_samples(sdb, f_pos_map)
        sdb.close()
        msg = '\n'.join(['barcodes:'] + ['%s -> %s' % (b, s.base_name()) for b, s in b2s.items()])
        self.logger.log(lg.DEBUG, msg)
        self.logger.log(lg.INFO, 'found %i samples.' % len(b2s))
        msg = 'features:\n' + '\n'.join('%s: %s' % (str(f), ','.join(str(x) for x in f.vals)) for f in features.values())
        self.logger.log(lg.DEBUG, msg)

        return bc_len, b2s, features

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
                if USER_STATES[cur] >= self.start_after:
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
        self.bw_dir = d + os.sep + self.bigwig_dirname
        self.bam_dir = d + os.sep + self.bam_dirname

        if os.path.isdir(self.tmp_dir): shutil.rmtree(self.tmp_dir)
        if os.path.islink(self.bw_dir): os.remove(self.bw_dir)
        elif os.path.isdir(self.bw_dir): shutil.rmtree(self.bw_dir)


        if self.start_after == BEGIN:
            # assuming all folder structure exists if check passes
            self.create_dir_and_log(d, lg.INFO)
            self.create_dir_and_log(self.fastq_dir)
            self.create_dir_and_log(self.bam_dir)
            self.create_dir_and_log(self.bw_dir)
            if self.keep_filtered:
                self.filtered_dir = d + os.sep + FILTERED_NAME
                self.create_dir_and_log(self.filtered_dir)
            if self.keep_unaligned:
                self.unaligned_dir = d + os.sep + UNALIGNED_NAME
                self.create_dir_and_log(self.unaligned_dir)

        if self.make_hub:
            if self.www_dir is None: self.www_dir = self.exp
            self.www_path = canonic_path('~/www/%s' % self.www_dir)

        self.create_dir_and_log(self.tmp_dir)

    def build_hub(self):
        self.logger.log(lg.INFO, 'Generating hub...')
        if os.path.isdir(self.www_path):
            self.logger.log(lg.DEBUG, 'Removing old folder')
            shutil.rmtree(self.www_path)
        self.create_dir_and_log(self.www_path)
        sacpath = self.www_path + os.path.sep + 'sacCer3'
        if not os.path.exists(sacpath): os.mkdir(sacpath)
        hubfile = open(self.www_path + os.path.sep + 'hub.txt', 'w')
        hubfile.write('\n'.join(["hub %s" % self.exp,
                                 "shortLabel %s" % self.exp,
                                 "longLabel %s" % self.exp,
                                 "genomesFile genomes.txt",
                                 "email %s" % self.hub_email]))
        genomesfile = open(self.www_path + os.path.sep + 'genomes.txt', 'w')
        genomesfile.write("genome sacCer3\n"
                          "trackDb sacCer3/trackDB.txt")
        trackfile = open(sacpath + os.path.sep + 'trackDB.txt', 'w')
        for s in self.samples.values():
            wurl = os.sep.join([URL_BASE+getpass.getuser(), self.www_dir,
                                self.bigwig_dirname, os.path.split(s.files['wbw'])[1]])
            curl = os.sep.join([URL_BASE+getpass.getuser(), self.www_dir,
                                self.bigwig_dirname, os.path.split(s.files['cbw'])[1]])
            hdr = '\n'.join(['track %s' % s.base_name(),
                             'container multiWig',
                             'aggregate transparentOverlay',
                             'type bigWig',
                             'autoScale on',
                             'visibility full',
                             'shortLabel %s' % s.base_name(),
                             'longLabel  %s' % s.full_name(),
                             'maxHeightPixels 100:32:8',
                             'priority 50'])
            wentry = '\n\t'.join(['\ttrack %s_w' % s.base_name(),
                                  'parent %s' % s.base_name(),
                                  'type bigWig 0 1000',
                                  'color 0,92,192',
                                  'alwaysZero on',
                                  'yLineOnOff on',
                                  'visibility full',
                                  'smoothingWindow 4',
                                  'windowingFunction mean',
                                  'graphTypeDefault bar',
                                  'bigDataUrl %s' % wurl])
            centry = '\n\t'.join(['\ttrack %s_c' % s.base_name(),
                                  'parent %s' % s.base_name(),
                                  'type bigWig 0 1000',
                                  'altcolor 0,92,192',
                                  'alwaysZero on',
                                  'yLineOnOff on',
                                  'visibility full',
                                  'smoothingWindow 4',
                                  'windowingFunction mean',
                                  'graphTypeDefault bar',
                                  'bigDataUrl %s' % curl])
            trackfile.write(hdr+'\n\n')
            trackfile.write(wentry+'\n\n')
            trackfile.write(centry+'\n\n')
        trackfile.close()
        mainurl = os.sep.join([URL_BASE+getpass.getuser(), self.www_dir, 'hub.txt'])
        new_path = self.www_path + os.sep + self.bigwig_dirname
        shutil.move(self.bw_dir, self.www_path)
        os.symlink(new_path, self.bw_dir, target_is_directory=True)
        sp.call('chmod -R 777 %s' % self.www_path, shell=True)
        msg = 'Transferred bigwig files to %s (link available in %s as well)' % (new_path, self.output_dir)
        self.logger.log(lg.DEBUG, msg)
        self.logger.log(lg.CRITICAL, 'Hub available at %s' % mainurl)

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
            cnt_path = self.tmp_dir + os.sep + (BC_COUNTS_FNAME % it)
            nobc = NO_BC_NAME + '-' + it
            arraydef = ';\n'.join('a["%s"]="%s-%s"' % (b, s, it) for b, s in b2s.items()) + ';\n'
            awk_str = (""" 'BEGIN {%s} {x=substr($4,1,%i); if (x in a) """,
                       """{c[a[x]]++; print >> "%s/"a[x];} else {c["%s"]++; print;} }""",
                       """END { for (bc in c) print bc, c[bc] >> "%s" } '""")
            awk_str = ''.join(awk_str) % (arraydef, self.bc_len , self.tmp_dir, nobc, cnt_path)
            return awk_str, cnt_path

        def merge_statistics(bc1, bc2):
            stat = "n_reads"
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
            for s in self.samples.values():
                if s.base_name() not in self.stats:
                    self.stats[s.base_name()][stat] += 0
            msg = '\n'.join(['%s: %i' % (s, c[stat]) for s, c in self.stats.items()])
            self.logger.log(lg.CRITICAL, 'read counts:\n' + msg)

        hb = {}
        for b,s in self.samples.items():
            hb.update({eb:s.base_name() for eb in hamming_ball(b, self.hamming_distance)})

        awk1p, cnt1 = compile_awk("1", {b: s.base_name() for b,s in self.samples.items()})
        awk2p, cnt2 = compile_awk("2", hb)
        outf = open(os.devnull, 'w') if no_bc is None else open(no_bc, 'wb')
        for r1, r2 in self.input_files:
            msg = 'splitting files:\n%s\n%s' % (os.path.split(r1)[1],os.path.split(r2)[1])
            self.logger.log(lg.INFO, msg)
            paste1 = sp.Popen('paste <(zcat %s) <(zcat %s)' % (r1,r2), stdout=sp.PIPE,
                              shell=True, executable='/bin/bash')
            awkin = sp.Popen(sh.split('paste - - - -'), stdin=paste1.stdout, stdout=sp.PIPE)
            if self.debug: # only a small subset of reads
                nlines = round(self.db_nlines/4)*4  # making sure it's in fastq units
                awkin = sp.Popen(sh.split('head -%i' % nlines), stdin=awkin.stdout, stdout=sp.PIPE)
            awk1 = sp.Popen(sh.split('awk -F "\\t" ' + awk1p), stdin=awkin.stdout, stdout=sp.PIPE)
            awk2 = sp.Popen(sh.split('awk -F "\\t" ' + awk2p), stdin=awk1.stdout, stdout=sp.PIPE)
            awkcmd = """awk -F "\\t" '{print $1"\\n"$3"\\n"$5"\\n"$7; print $2"\\n"4"\\n"$6"\\n"$8;}' """
            wfastq = sp.Popen(sh.split(awkcmd), stdin=awk2.stdout, stdout=sp.PIPE)
            gzip = sp.Popen(['gzip'], stdin=wfastq.stdout, stdout=outf)
            wfastq.wait() # to prevent data interleaving
        gzip.wait()
        self.logger.log(lg.INFO, 'Barcode splitting finished.')

        merge_statistics(cnt1, cnt2)

    def aftermath(self):
        # remove temp folder
        # modify file permissions for the entire tree
        # make everything read only (optional?)
        # store pipeline code?
        # merge and report statistics
        # merge data to single (usable) files
        if self.debug is None:
            shutil.rmtree(self.tmp_dir)

        pickle.dump(self.args, open(self.output_dir + os.sep + 'args.pkl', 'wb'))

        # change permissions so everyone can read into folder
        for d, _, fs in os.walk(self.output_dir):
            st = os.stat(d)
            os.chmod(d, st.st_mode | stat.S_IRGRP | stat.S_IXGRP)

        self.logger.log(lg.CRITICAL, 'All done.')
        self.copy_log()

    def create_dir_and_log(self, path, level=lg.DEBUG):
        create_dir(path)
        self.logger.log(level, 'creating folder %s' % path)

    def add_hub(self):
        pass
        # track_db = open(
        # for (hf, hr), p, hd in zip(handles, paths, hdr):
        #     hd + '_F', URL + hd + '_F.bw', 'pA - %s (F)' % hd, 'Wilkening 2013 polyA data, %s (forward)' % hd)
        #     trackfile.write("track %s\n"
        #                     "bigDataUrl %s\n"
        #                     "shortLabel %s\n"
        #                     "longLabel %s\n"
        #                     "type bigWig\n"
        #                     "visibility full\n"
        #                     "viewLimits 0:500\n\n" % entries)
        #     entries = (
        #     hd + '_R', URL + hd + '_R.bw', 'pA - %s (R)' % hd, 'Wilkening 2013 polyA data, %s (reverse)' % hd)
        #     trackfile.write("track %s\n"
        #                     "bigDataUrl %s\n"
        #                     "shortLabel %s\n"
        #                     "longLabel %s\n"
        #                     "type bigWig\n"
        #                     "visibility full\n"
        #                     "viewLimits -500:0\n\n" % entries)
        # hubfile = open(hub_path + os.path.sep + 'hub.txt', 'wb')
        # hubfile.write("hub %s\n"
        #               "shortLabel pA sites\n"
        #               "longLabel Data relevant to 3' processing and alternative UTRs\n"
        #               "genomesFile genomes.txt\n"
        #               "email  alonappleboim@gmail.com" % os.path.split(hub_path)[-1])
        # genomesfile = open(hub_path + os.path.sep + 'genomes.txt', 'wb')
        # genomesfile.write("genome sacCer3\n"
        #                   "trackDb sacCer3%strackDB.txt" % os.path.sep)


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
                   help='path to the folder in which most files are written. '
                        'If not given, the date and info from barcode file are used to '
                        'generate a new folder in %s' % DATA_PATH)
    g.add_argument('--fastq_dirname', '-fd', default='FASTQ', type=str,
                   help='name of folder in which fastq files are written. relative to output_dir')
    g.add_argument('--bam_dirname', '-bd', default='BAM', type=str,
                   help='name of folder in which bam files are written relative to output_dir')
    g.add_argument('--bigwig_dirname', '-bwd', default='BIGWIG', type=str,
                   help='name of folder in which bigwig files are written relative to output_dir')
    g.add_argument('--user_emails', '-ue', default=None, type=str,
                   help="if provided these comma separated emails will receive notifications of ctitical "
                        "events (checkpoints, fatal errors, etc.)")
    g.add_argument('--debug', '-d', default=None, type=str,
                   help='Highly recommended. Use this mode with a pair of comma separated integers:'
                        '<numlines>,<numsamples>. The pipeline will extract this number of lines from '
                        'every input file pair, and will only run with this number of samples out of the '
                        'given barcode file. If output folder (-od) is not given, results are written to '
                        'a "debug" folder.')

    g = p.add_argument_group('Barcode Splitting')
    g.add_argument('--sample_db', '-sd', type=str, default=None,
                   help='a file with an experiment name in first row \n experiment: <expname> followed by a'
                        ' header whose first column is "barcode", and the remaining entries are the features '
                        'present in the experiment: <name>(short_name):(str|int|float)[units], with the '
                        'units and short_name being optional. The remaining lines define the experiment samples '
                        'according to the header. Default is "sample_db.csv" in the folder of the --fastq_prefix input')
    g.add_argument('--umi_length', '-ul', type=int, default=8,
                   help='UMI length')
    g.add_argument('--hamming_distance', '-hd', default=1, type=int,
                   help='barcode upto this hamming distance from given barcodes are handled by '
                        'the pipeline')
    g.add_argument('--keep_nobarcode', '-knb', action='store_true',
                   help='keep reads that did not match any barcode in a fastq file. Notice that the '
                        'R1/R2 reads are interleaved in this file.')

    g = p.add_argument_group('Alignment')
    g.add_argument('--scer_index_path', '--sip', type=str, default='/cs/wetlab/genomics/scer/bowtie/sacCer3',
                   help='path prefix of s. cervisae genome bowtie index')
    g.add_argument('--klac_index_path', '-kip', type=str, default=None,
                   help='If given, data is also aligned to k. lacis genome (can be found in '
                        '/cs/wetlab/genomics/klac/bowtie/genome)')
    g.add_argument('--n_threads', '-an', type=int, default=4,
                   help='number of threads used for alignment per bowtie instance')
    g.add_argument('--keep_unaligned', '-ku', action='store_true',
                   help='if set, unaligned reads are written to '
                        'output_folder/%s/<sample_name>.bam' % UNALIGNED_NAME)

    g = p.add_argument_group('Filter',
                             description='different filters applied to base BAM file. Only reads that pass all filters '
                                         'are passed on')
    g.add_argument('--filter', '-F', action='store',
                   default='dup(),qual()',
                   help='specify a filter scheme to apply to data. Expected string conforms to:\n' \
                        '[<filter_name>([<argname1=argval1>,]+)[+|-]);]*\n. use "run -fh" for more info. '
                        'default = "dup(),qual()"')
    g.add_argument('--keep_filtered', '-kf', action='store_true',
                   help='if set, filtered reads are written to '
                        'output_folder/%s/<sample_name>.bam' % FILTERED_NAME)
    g.add_argument('--filter_specs', '-fh', action='store_true',
                   help='print available filters, filter help and exit')

    g = p.add_argument_group('Execution and third party')
    g.add_argument('--exec_on', '-eo', choices=['slurm', 'local'], default='local', #TODO implement slurm...
                   help='whether to to submit tasks to a slurm cluster (default), or just use the local processors')
    g.add_argument('--bowtie_exec', '-bx', type=str, default=BOWTIE_EXEC,
                   help='full path to bowtie executable')
    g.add_argument('--samtools_exec', '-sx', type=str, default=SAMTOOLS_EXEC,
                   help='full path to samtools executable')

    g = p.add_argument_group('Tracks')
    g.add_argument('--make_hub', '-mh', action='store_false',
                   help='prevent the pipeline from genrating a browser hub in your www folder')
    g.add_argument('--hub_email', '-he', action='store', default='noemail@nodomain.com',
                   help='the contact email for the generated hub')
    g.add_argument('--www_dir', '-wd', action='store', default=None,
                   help='the directory in which the hub is generated. If not given, '
                        'experiment name is used.')

    g = p.add_argument_group('Count')
    g.add_argument('--tts_file', '-tf', default=TTS_MAP,
                   help='annotations for counting. Expected format is a tab delimited file with "chr", "ACC", "start",'
                        '"end", and "TTS" columns. default is found at %s' % TTS_MAP)
    g.add_argument('--keep_tts_bed', '-ktb', action='store_true',
                   help='whether the TTS window definition bed file should be kept or dicarded')
    g.add_argument('--dont_bound_start', '-dbs', action='store_false',
                   help='do not limit counting window orf start (i.e. if longer thatn ORF, trim counting window)')
    g.add_argument('--count_window', '-cw', type=str, default='[-750,250]',
                   help='Comma separated limits for tts counting, relative to annotation TTS. default=-100,500')

    g = p.add_argument_group('Export')
    g.add_argument('--exporters', '-E', action='store',
                   default='tab();mat(r=True)',
                   help='specify a exporters for the pipeline data and statistics. default = "tab();mat(r=True)"')
    g.add_argument('--export_path', '-ep', default=None,
                   help='if given, exported data is copied to this path as well')
    g.add_argument('--exporter_specs', '-eh', action='store_true',
                   help='print available exporters, exporter help and exit')
    return p


def pprint_class_with_args_dict(dict):
    slist = []
    for cls in dict.values():
        args = cls.__dict__['args']
        slist.append('Name: %s' % cls.__dict__['name'])
        for a, v in cls.__dict__.items():
            if a.startswith('__'): continue
            if hasattr(v, '__call__'): continue
            if a not in ['args', 'name']:
                slist.append('\t%s: %s' % (str(a),str(v)))
        slist.append('\targuments:')
        for name, (type, default, desc) in args.items():
            t = 'int' if type is int else 'float' if type is float else 'bool' if type is bool else 'str'
            slist.append('\t\t%s(%s):, default: %s, description: %s' % (name, t, str(default), desc))
    return '\n'.join(slist)


def parse_args(p):
    """
    :param p: argument parser
    :return: the arguments parsed, after applying all argument logic and conversion
    """
    args = p.parse_args()

    if args.filter_specs:
        h = ('Any collection of filters can be applied. The only reads that are written to the BAM '
             'folder are the ones that pass the complete list of filters. If you want to keep filtered '
             'reads for inspection, use the --keep_filtered options, and they will be written to the '
             '%s folder. Note that the reads are written by each filter separately, so any read that was '
             'filtered by multiple filters will be written more than once. \n'
             'A collection of filters is given as a comma seprated list of filter calls - a filter '
             'name followed by parentheses with optional argument setting within. The parentheses are followed '
             'by an optional +|- sign, to negate the filter ("-"). The default is "+". For example:\n'
             '"dup(kind=start&umi),polya(n=5)" will only keep reads that are unique when considering only read '
             'start position and the umi, and are considered polya reads with more than 5 A/T. See all available '
             'filters below.\n') % (FILTERED_NAME,)
        spec = pprint_class_with_args_dict(collect_filters())
        print('\n========= FILTER HELP and SPECIFICATIONS =========\n')
        print(h)
        print(spec)
        exit()

    if args.exporter_specs:
        h = ('Any collection of exporters can be given. Exporters format the quantiative output of the piplines '
             'for downstream analysis. Specifically - the statistics and TTS read counts. In addition,'
             'using the --export_path option will copy the outputs to the requested folder with a prefix of the '
             'experiment name.\n'
             'A collection of exporters is given as a semicolon-separated list of exporter calls - an exporter name '
             'followed by parentheses with optional argument setting within. For example:\n'
             '"tab(),mat(r=True)" export all files in a tab delimited format and a matlab format with table reshaping'
             ' according to experimental features.')
        spec = pprint_class_with_args_dict(collect_exporters())
        print('\n========= EXPORTER HELP and SPECIFICATIONS =========\n')
        print(h)
        print(spec)
        exit()

    args.__dict__['start_after'] = USER_STATES[args.start_after]
    if args.start_after != BEGIN:
        if args.output_dir is None:
            print('If the --start_after option is used, an existing output directory must be provided (-od).')
            exit()
        args.__dict__['fastq_prefix'] = args.output_dir  # ignoring input folder
    else:
        if args.fastq_prefix is None:
            print('If the --start_after option is not used, an input '
                  'fastq prefix/folder must be provided (--fastq_prefix).')
            exit()
    p, s = os.path.split(args.fastq_prefix)
    args.__dict__['fastq_path'] = p
    args.__dict__['fastq_pref'] = s
    args.__dict__['fastq_prefix'] = canonic_path(args.fastq_prefix)
    p, s = os.path.split(args.fastq_prefix)

    if args.sample_db is None:
        args.__dict__['sample_db'] = os.sep.join([args.fastq_path, 'sample_db.csv'])

    if args.debug is not None:
        nlines, nsamples = args.debug.split(',')
        args.__dict__['db_nlines'] = int(nlines)
        args.__dict__['db_nsamples'] = int(nsamples)

    if args.output_dir is not None:
        args.__dict__['output_dir'] = canonic_path(args.__dict__['output_dir'])

    if args.export_path is not None:
        args.__dict__['export_path'] = canonic_path(args.export_path)

    args.__dict__['count_window'] = [int(x) for x in args.count_window[1:-1].split(',')]
    return args


if __name__ == '__main__':
    p = build_parser()
    a = parse_args(p)
    mh = MainHandler(a, ' '.join(sys.argv)) #also print git current version from ./.git/log/HEAD,last row, 2nd column
    mh.execute()