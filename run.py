"""
This script performs the preprocessing steps and starts a dedicated process for every sample.
 """


import threading
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
import traceback
from queue import Empty

from config import *

STATES = {
    0:'BEGINNING',
    1:'FASTQ',
    2:'ALIGN',
    3:'FILTER',
}
STATE_ORDER = {
    'BEGINNING':0,
    'FASTQ':1,
    'ALIGN':2,
    'FILTER':3,
}

# global functions should not log to logger, wrap them with
# instance functions if you want them to log things (only) through logq

def create_dir(path):
    if not os.path.exists(path): os.mkdir(path)


def get_logfile():
    now = datetime.datetime.now()
    p = LOG_PATH + os.sep + str(now.year)
    create_dir(p)
    p += os.sep + str(now.month)
    create_dir(p)
    fname = '%s_%i_%i-%i-%i' % (getpass.getuser(), now.day, now.hour, now.minute, now.second)
    return p + os.sep + fname


def parse_bowtie_stats(bt_stats):
    # parse this output:
    # <line>*
    # 10263 reads;
    # of these:
    #     10263(100.00 %) were unpaired;
    #     of these:
    #         1055(10.28 %) aligned 0 times
    #         4337(42.26 %) aligned exactly 1 time
    #         4871(47.46 %) aligned > 1 times
    # 89.72 % overall alignment rate
    # <line>*
    nlines, Ns = 0, []
    for i, line in enumerate(bt_stats):
        m = re.match('\s*(\d+).*', line)
        if m is not None:
            nlines += 1
            if nlines in [2,3,4,5]:
                Ns.append(m.group(1))
    return Ns


def hamming_ball(seq, radius, alphabet='CGTAN'):
    ball = [seq]
    if radius > 0:
        for i in range(len(seq)):
            for l in alphabet:
                seql = list(seq)
                seql[i] = l
                ball.extend(hamming_ball(''.join(seql), radius-1, alphabet))
    return list(set(ball))


def isint(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def block_exec(command, on_slurm, stdout=None):
    if on_slurm:
        p = sp.Popen(['srun', command], stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = p.communicate()
    else:
        p = sp.Popen(sh.split(command))
        out, err = p.communicate()
    return out, err


def canonic_path(fname):
    return os.path.abspath(os.path.expanduser(fname))


def write_logs(logq, user_emails=None):

    #configure logging behavior
    logger = lg.getLogger()
    logfile = get_logfile()
    logger.setLevel(lg.DEBUG)
    fh = lg.FileHandler(logfile)  # copy log to the log folder when done.
    fh.formatter = lg.Formatter('%(levelname)s\t%(asctime)s\t%(message)s')
    ch = lg.StreamHandler()
    ch.formatter = lg.Formatter('%(asctime)-15s\t%(message)s')
    fh.setLevel(lg.DEBUG)
    ch.setLevel(lg.INFO)
    if user_emails is not None:
        mailh = ThreadedTlsSMTPHandler(mailhost=('smtp.gmail.com', 587),
                                       fromaddr='transeq.pipeline@google.com',
                                       toaddrs=user_emails.split(','),
                                       credentials=('transeq.pipeline', 'transeq1234'),
                                       subject='TRANSEQ pipeline message')
        mailh.setLevel(lg.CRITICAL)
        logger.addHandler(mailh)
    logger.addHandler(fh)
    logger.addHandler(ch)

    # get logs from Q and write them...
    for lvl, msg in iter(logq.get, None):
        logger.log(lvl, msg)

    # send logfile name
    logq.put(logfile)


class MainHandler(object):

    def __init__(self, argobj, cmdline):
        self.__dict__.update(argobj.__dict__)

        self.start_from = self.start_after

        self.fastqq = mp.Queue()
        self.bamq = mp.Queue()
        self.filterq = mp.Queue()
        self.analysisq = mp.Queue()

        self.logq = mp.Queue()
        self.logwrtr = mp.Process(target=write_logs, args=(self.logq, self.user_emails))
        if self.debug:
            self.logq.put((lg.INFO, '=================== DEBUG MODE (%s)===================' % self.debug))
        self.logwrtr.start()

        self.check_third_party()
        self.parse_filters()
        self.parse_analyses()

    def create_dir_and_log(self, path, level=lg.DEBUG):
        create_dir(path)
        self.logq.put((level, 'creating folder %s' % path))

    def start_worker(self, cls, args):
        worker = cls(self.__dict__, *args)
        worker.daemon = True
        worker.start()
        return worker

    def execute(self):
        self.parse_barcode_file()
        self.generate_folder_structure()
        self.write_barcode_file()

        fastq_ncomplete = set([])
        if self.start_from == STATES[0]:
            self.collect_input_fastqs()
            self.split_barcodes()
            for bc, sample in self.b2s.items():
                fastq_ncomplete.add(sample)
                self.start_worker(FormatFastq, (sample, self.logq, self.fastqq))

        bam_ncomplete = set([])
        while fastq_ncomplete:
            sample = self.fastqq.get()
            fastq_ncomplete.remove(sample)
            bam_ncomplete.add(sample)
            self.start_worker(Fastq2BAM, (sample, self.logq, self.bamq))
        if self.start_from == STATES[0]: self.mark_fastq()

        if self.start_from == STATES[1]:
            for bc, sample in self.b2s.items():
                bam_ncomplete.add(sample)
                self.start_worker(Fastq2BAM, (sample, self.logq, self.bamq))
        else:
            assert bam_ncomplete

        filter_ncomplete = set([])
        while bam_ncomplete:
            sample = self.bamq.get()
            bam_ncomplete.remove(sample)
            for fname, filter in self.filters.items():
                filter_ncomplete.add((sample, fname))
                self.start_worker(FilterBAM, (sample, filter, self.logq, self.filterq))
        if STATE_ORDER[self.start_from] <= 1: self.mark_aligned()

        # if self.start_from == STATES[3]:
        #     for bc, sample in self.b2s.items():
        #         for fname, filter in self.filters.items():
        #             filter_ncomplete.add((sample, fname))
        #             self.start_worker(FilterBAM, (sample, filter, self.logq, self.filterq))
        #
        # analysis_ncomplete = set([])
        # while filter_ncomplete:
        #     sample, fname = self.filterq.get()
        #     filter_ncomplete.remove((sample, fname))
        #     for aname, analysis in self.anlyses.items():
        #         analysis_ncomplete.add((sample, fname, aname))
        #         self.start_worker(AnalyzeBAM, (sample, fname, analysis, self.logq, self.analysisq))
        # if self.start_from == STATES[3]: self.mark_filtered()
        #
        # if self.start_from == STATES[4]:
        #     for bc, sample in self.b2s.items():
        #         for fname, filter in self.filters.items():
        #             for aname, analysis in self.analyses.items():
        #                 analysis_ncomplete.add((sample, fname, aname))
        #                 self.start_worker(AnalyzeBAM, (sample, fname, analysis, self.logq, self.analysisq))
        #
        # while analysis_ncomplete:
        #     sample, fname, aname = self.analysisq.get()
        #     analysis_ncomplete.remove((sample, fname, aname))
        #
        # # finally, merge results
        # mergers = []
        # for aname, analysis in self.analyses.items():
        #     mergers.append(self.start_worker(MergeAnalysis, (aname, analysis, self.logq)))
        # for m in mergers: m.join()

        # clean up and go home
        self.aftermath()

    def print_stats(self):
        stats = set([])
        for s, st in self.stats.items():
            stats |= set(st.keys())
        fns = ['sample'] + list(stats)
        fh = open(self.output_dir + os.sep + 'sample_statistics.csv','w')
        wrtr = csv.DictWriter(fh, fns)
        wrtr.writeheader()
        for s, st in self.stats.items():
            wrtr.writerow(dict(sample=s, **{k:str(st[k]) for k in fns[1:]}))
        fh.close()

    def mark(self, stage):
        self.start_after = stage
        self.print_stats()
        msg = ('Finished stage %s. You can continue the pipeline from this point '
               'with the option -sa %s (--start_after %s)'  % (stage, stage, stage))
        self.logq.put((lg.INFO, msg))
        fh = open(self.output_dir + os.sep + '.pipeline_state', 'w')
        fh.write('%s\n' % stage)
        fh.close()

    def get_mark(self):
        fh = open(self.output_dir + os.sep + '.pipeline_state')
        stage = fh.read().strip()
        fh.close()
        return stage

    def mark_split(self):
        """
        collect all statistics upto this point and mark folder with apropriate flags
        for future easy reference (see "start_after" option of this script)
        """
        pass

    def mark_fastq(self):
        """
        collect all statistics upto this point and mark folder with apropriate flags
        for future easy reference (see "start_after" option of this script)
        """
        self.logq.put((lg.INFO, 'All fastq files were written to: %s' % self.fastq_dir))
        self.mark(STATES[1])

    def mark_aligned(self):
        # collect all statistics upto this point and mark folder with apropriate flags
        # for future easy reference (see "start_after" option of this script)
        for f in os.listdir(self.tmp_dir):
            fs = f.split(os.extsep)
            if fs[-1]!= 'stats': continue
            if fs[1] == 'bowtie':
                sample = fs[0]
                f = self.tmp_dir + os.sep + f
                with open(f) as F:
                    for line in F:
                        stat, cnt = line.strip().split('\t')
                        if len(fs) == 4: stat += '_'+fs[2]
                        self.stats[sample][stat] += int(cnt)
                os.remove(f)
        self.mark(STATES[2])

    def mark_filtered(self):
        # collect all statistics upto this point and mark folder with apropriate flags
        # for future easy reference (see "start_after" option of this script)
        self.mark(STATES[4])
        pass

    def mark_analyzed(self):
        # collect all statistics upto this point and mark folder with apropriate flags
        # for future easy reference (see "start_after" option of this script)
        pass

    def prepare_outputs(self):
        pass

    def check_third_party(self):
        if self.exec_on == 'slurm':
            try:
                p = sp.Popen(sh.split('srun "slurm, are you there?"'), stdout=sp.PIPE, stderr=sp.PIPE)
                p.communicate()
                self.logq.put((lg.INFO, "slurm check.. OK"))
            except OSError as e:
                self.logq.put((lg.CRITICAL, "This is not a slurm cluster, execute with flag -eo=local"))
                raise e

        for ex in [k for k in self.__dict__.keys() if k.endswith('exec')]:
            try:
                p = sp.Popen(sh.split('%s --help' % self.__dict__[ex]), stdout=sp.PIPE, stderr=sp.PIPE)
                p.communicate()
                self.logq.put((lg.INFO, "%s check.. OK" % ex))
            except OSError as e:
                self.logq.put((lg.CRITICAL, "could not resolve %s path: %s" % (ex, self[ex])))
                raise e

    def parse_barcode_file(self):
        with open(self.barcode_file) as BCF:
            name = re.match('.*name.*:\s+(\w+)', BCF.readline())
            exp = re.match('.*experiment.*:\s+(\w+)', BCF.readline())
            if name is None or exp is None:
                msg = 'barcodes file should contain a 2 line header with name: ' \
                      '<name> and experiment: <expname> (no spaces)'
                self.logq.put((lg.CRITICAL, msg))
                raise ValueError(msg)
            self.name = name.group(1)
            self.exp = exp.group(1)
            msg = 'name: %s, experiment: %s' % (self.name, self.exp)
            self.logq.put((lg.INFO, msg))
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
            self.logq.put((lg.DEBUG, msg))
            self.logq.put((lg.INFO, 'found %i samples.' % len(b2s)))

        self.b2s = b2s
        self.s2b = s2b
        self.stats = {s: Counter() for s in s2b.keys()}

    def write_barcode_file(self):
        shutil.copy(self.barcode_file, self.output_dir + os.sep + 'barcodes')

    def collect_input_fastqs(self):
        files = {}

        for fn in os.listdir(self.fastq_path):
            if not fn.endswith('fastq.gz'): continue
            if not fn.startswith(args.fastq_pref): continue
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
        self.logq.put((lg.INFO, msg))
        if not files:
            msg = "could not find R1/R2 fastq.gz pairs in given folder: %s" % self.fastq_path
            self.logq.put((lg.CRITICAL, msg))
            raise IOError(msg)
        self.input_files = files

    def generate_folder_structure(self):

        if self.start_from != STATES[0]:
            try:
                cur = self.get_mark()
                if STATE_ORDER[cur] >= STATE_ORDER[self.start_after]:
                    msg = 'restarting from %s in folder: %s ' % (self.start_after, self.output_dir)
                    self.logq.put((lg.INFO, msg))
                else:
                    msg = 'folder state %s in folder %s incompatible with start_after %s request' \
                          % (cur, self.output_dir, self.start_after)
                    self.logq.put((lg.CRITICAL, msg))
                    exit()
            except IOError:
                msg = 'could not find an existing output folder: %s' % self.output_dir
                self.logq.put((lg.CRITICAL, msg))
                exit()
        else:
            if self.output_dir is None:
                folder = canonic_path(DATA_PATH) + os.sep + self.name[0].upper() + self.name[1:].lower()
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

        if self.start_from == STATES[0]:
            # assuming all folder structure exists if check passes
            self.create_dir_and_log(d, lg.INFO)
            self.create_dir_and_log(self.fastq_dir)
            self.create_dir_and_log(self.bam_dir)
            self.create_dir_and_log(self.tmp_dir)
            if self.keep_unaligned:
                self.create_dir_and_log(d + os.sep + UNALIGNED_NAME)

    def split_barcodes(self):
        #
        #  first - compile awk scripts that will split input according to barcodes,
        # then use them with some smart pasting. For each input pair:
        # paste <(zcat {R1}) <(zcat {R2}) |\ %paste both files one next to the other
        # paste - - - -|\ % make every 4 lines into one
        # awk -F "\\t" -f {script1-exact_match}' |\ % split exact barcode matches
        # awk -F "\\t" -f {script2-ham_dist}' % split erroneous matches (according to hamming distance)
        #

        def compile_awk(it, b2s):
            awk_path = self.tmp_dir + os.sep + 'split_barcodes_%s.awk' % it
            cnt_path = self.tmp_dir + os.sep + 'bc_counts-%s' % it
            nobc = 'no-barcode-' + it
            script = open(awk_path, 'w')
            arraydef = ';\n'.join('a["%s"]="%s-%s"' % (b, s, it) for b, s in b2s.items()) + ';\n'
            script.write("""
                BEGIN {
                %s
                }
                {
                x=substr($4,1,7); if (x in a) {c[a[x]]++; print >> "%s/"a[x];} else {c["%s"]++; print;}
                }
                END {
                for (bc in c)
                print bc, c[bc] >> "%s"
                }
                """ % (arraydef, self.tmp_dir, nobc, cnt_path))
            script.close()
            return awk_path

        def merge_statistics():
            bc1 = os.sep.join([self.tmp_dir, "bc_counts-1"])
            bc2 = os.sep.join([self.tmp_dir, "bc_counts-2"])
            stat = "#reads"
            self.stats['no-barcode'] = Counter()
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
                    no_counts.append(s)
                    self.stats[s][stat] += 0
            msg = '\n'.join(['%s: %i' % (s, c['#reads']) for s, c in self.stats.items()])
            self.logq.put((lg.CRITICAL, 'read counts:\n' + msg))

        hb = {}
        for b,s in self.b2s.items():
            hb.update({eb:s for eb in hamming_ball(b, self.hamming_distance)})

        awk1p = compile_awk("1", self.b2s)
        awk2p = compile_awk("2", hb)
        for r1, r2 in self.input_files:
            msg = 'splitting files:\n%s\n%s' % (os.path.split(r1)[1],os.path.split(r2)[1])
            self.logq.put((lg.INFO, msg))
            paste1 = sp.Popen('paste <(zcat %s) <(zcat %s)' % (r1,r2), stdout=sp.PIPE,
                              shell=True, executable='/bin/bash')
            awkin = sp.Popen(sh.split('paste - - - -'), stdin=paste1.stdout, stdout=sp.PIPE)
            if self.debug: # only a small subset of reads
                awkin = sp.Popen(sh.split('head -%i' % self.db_nlines), stdin=awkin.stdout, stdout=sp.PIPE)
            awk1 = sp.Popen(sh.split('awk -F "\\t" -f ' + awk1p), stdin=awkin.stdout, stdout=sp.PIPE)
            awk2 = sp.Popen(sh.split('awk -F "\\t" -f ' + awk2p), stdin=awk1.stdout, stdout=open(os.devnull,'w'))
            awk2.wait()
            self.logq.put((lg.INFO, 'Barcode splitting finished.'))

    def parse_filters(self):
        self.filters = {}
        # try:
        #     self.filters.split(';')
        # except AttributeError as e:
        #     self.logq.put((lg.ERROR), 'could not resolve filter argument, see -fh')
        #     raise e

    def parse_analyses(self):
        self.analyses = {}

    def aftermath(self):
        # remove temp folder
        # make everything read only (optional?)
        # store pipeline code?
        # merge and report statistics
        # merge data to single (usable) files

        self.logq.put((lg.CRITICAL, 'All done.'))
        self.logq.put_nowait(None)
        self.logwrtr.join()
        logfile = self.logq.get()
        shutil.copy(logfile, self.output_dir + os.sep + 'full.log')


class FormatFastq(mp.Process):

    def __init__(self, argdict, sample, logq, comq):
        mp.Process.__init__(self)
        self.sample = sample
        self.logq = logq
        self.comq = comq
        self.__dict__.update(argdict)

    def build_fastq(self):
        # compiling a simple bash script to compile a minimal fastq.gz file:
        #
        # cat {split}-1 {split}-2 |\
        # awk -F "\\t" '{print "@umi:"substr($4,8,4)"\\n"$3"\\n+\\n"$7}'|\
        # gzip > {fastq.gz}
        #
        fname = self.sample + '.fastq.gz'
        target_path = self.fastq_dir + os.sep + fname
        pre1, pre2 = self.tmp_dir+os.sep+self.sample+'-1', self.tmp_dir+os.sep+self.sample+'-2'
        if os.path.isfile(pre1) and os.path.isfile(pre1):
            cat = sp.Popen(['cat', pre1, pre2], stdout=sp.PIPE)
        elif os.path.isfile(pre1):
            cat = sp.Popen(['cat', pre1], stdout=sp.PIPE)
        elif os.path.isfile(pre2):
            cat = sp.Popen(['cat', pre2], stdout=sp.PIPE)
        else:
            cat = sp.Popen(['cat'], stdin=os.devnull, stdout=sp.PIPE)
        awk = sp.Popen(sh.split('''awk -F "\\t" '{print "@umi:"substr($4,8,4)"\\n"$3"\\n+\\n"$7}' '''),
                       stdin=cat.stdout, stdout=sp.PIPE)
        gzip = sp.Popen(['gzip'], stdin=awk.stdout, stdout=open(target_path,'wb'))
        gzip.wait()
        if os.path.isfile(pre1): os.remove(pre1)
        if os.path.isfile(pre2): os.remove(pre2)
        self.logq.put((lg.DEBUG, '%s. ready' % fname))

    def run(self):
        try:
            self.logq.put((lg.DEBUG, 'Started working on %s' % self.sample))
            self.build_fastq()
            self.comq.put(self.sample) #fastq ready

        except Exception as e:
            self.logq.put((lg.CRITICAL, 'process ended abruptly:\n%s' % traceback.format_exc(10)))
            raise e


class Fastq2BAM(mp.Process):

    def __init__(self, argdict, sample, logq, comq):
        mp.Process.__init__(self)
        self.sample = sample
        self.logq = logq
        self.comq = comq
        self.__dict__.update(argdict)
        self.fastq = self.fastq_dir + os.sep + self.sample + '.fastq.gz'

    def fastq2bam(self):
        # align, bam, sort and index:
        #
        # bowtie2 --local -p {N} -U {fastq.gz} -x {index} 2> {stats} | samtools view -b > {tmp}
        # samtools view -F4 -b | samtools sort > {bam}
        # samtools index {bam}
        # rm {tmp}
        #
        #
        tmpbam = self.tmp_dir + os.sep + self.sample + '.tmp.bam'
        tmpstats = self.tmp_dir + os.sep + self.sample + '.bowtie.stats'
        fname = self.sample + '.bam'
        self.bam = os.sep.join([self.bam_dir, fname])
        nbam = os.sep.join([self.output_dir, UNALIGNED_NAME, fname])
        bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' %
                               (self.bowtie_exec, self.n_threads, self.fastq, self.scer_index_path)),
                      stdout=sp.PIPE, stderr=sp.PIPE)
        st = sp.Popen(sh.split('samtools view -b -o %s' % tmpbam), stdin=bt.stdout)
        st.wait()
        if self.keep_unaligned:
            naligned = sp.Popen(sh.split('samtools view -f4 -b %s -o %s' % (tmpbam, nbam)), stdout=sp.PIPE)
            self.logq.put((lg.INFO, 'unaligned BAM saved to: %s' % nbam))
            naligned.wait()
        aligned = sp.Popen(sh.split('samtools view -F4 -b %s' % tmpbam), stdout=sp.PIPE)
        sort = sp.Popen(sh.split('samtools sort -o %s' % self.bam), stdin=aligned.stdout)
        sort.wait()
        ind = sp.Popen(sh.split('samtools index %s' % self.bam))
        ind.wait()
        os.remove(tmpbam)
        Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
        msg = '%s aligned. Total: %s, No: %s, Unique: %s, Non-Unique: %s)' % ((self.sample,) + tuple(Ns))
        self.logq.put((lg.INFO, msg))
        with open(tmpstats, 'w') as S:
            S.write('\n'.join(['total\t%s' % Ns[0],
                               'no-align\t%s' % Ns[1],
                               'unique-align\t%s' % Ns[2],
                               'multiple-align\t%s' % Ns[3]]))

    def klac_count(self):
        # align, and parse statisticsbam, sort and index.
        # bowtie2 --local -p 4 -U {fastq.gz} -x {index} 2> {stats} >/dev/null
        if self.klac_index_path is None: return
        tmpstats = self.tmp_dir + os.sep + self.sample + '.bowtie.klac.stats'
        bt = sp.Popen(sh.split('%s --local -p 4 -U %s -x %s' % (self.bowtie_exec, self.fastq, self.klac_index_path)),
                      stderr=sp.PIPE, stdout=open(os.devnull, 'w'))
        Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
        msg = '%s aligned to k. lcatis. Total: %s, No: %s, Unique: %s, Non-Unique: %s)' % ((self.sample,) + tuple(Ns))
        self.logq.put((lg.INFO, msg))
        with open(tmpstats, 'w') as S:
            S.write('\n'.join(['total\t%s' % Ns[0],
                               'no-align\t%s' % Ns[1],
                               'unique-align\t%s' % Ns[2],
                               'multiple-align\t%s' % Ns[3]]))

    def run(self):
        try:
            ct = threading.Thread(target=self.klac_count)
            ct.start()
            self.fastq2bam()
            ct.join()  # making sure lactis alignment is also complete
            self.comq.put(self.sample) #bam complete

        except Exception as e:
            self.logq.put((lg.CRITICAL, 'process ended abruptly:\n%s' % traceback.format_exc(10)))
            raise e


class FilterBAM(mp.Process):
    pass


class AnalyzeBAM(mp.Process):
    pass


class MergeAnalysis(mp.Process):
    pass


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('--fastq_prefix', '-fp', type=str, default=None,
                   help='path to a prefix of fastq files (R1 & R2) containing the transeq data.'
                        'This can be a folder (must end with "/"), in which case all R1/R2 pairs'
                        'in the folder are considered, or a "path/to/files/prefix", in which case '
                        'all files in the path with the prefix are considered')
    g.add_argument('--start_after', '-sa', default=STATES[0],
                   choices=[STATES[k] for k in sorted(STATE_ORDER.values())],
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
    g.add_argument('--hamming_distance', '-hd', default=1, type=int,
                   help='barcode upto this hamming distance from given barcodes are handled by '
                        'the pipeline')

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

    # filters = get_filters()
    g = p.add_argument_group('Filters',
                             description='different filters applied to base BAM file. Each filter result is '
                                         'processeed downstream and reported separately')
    g.add_argument('--filters', '-F', action='store',
                   default='unique:dup()+,qual(maxq=1,)+;unique-polya:dup()+,qual(1,)+,polya(5,)+',
                   help='specify filter schemes to apply to data. Expected string conforms to ([] are for grouping):\n' \
                        '[<filter_scheme>:<filter>([<argname1=argval1>,]+)[+|-]);]*\n the filter_scheme will be used to name all '
                        'resulting outputs from this branch of the data. use "run -fh" for more info.')
    g.add_argument('--filter_specs', '-fh', action='store_true',
                   default='print available filter specs and filter help and exit')

    g = p.add_argument_group('Analyses',
                             description='different analyses applied to filtered BAM files. Each analysis is '
                                         'performed on single samples, and when all is done, results are merged '
                                         'to a single file in the output/path/results folder')
    g.add_argument('--analyses', '-A', action='store',
                   default='coverage(); XXXX',
                   help='XXX')
    g.add_argument('--analysis_specs', '-ah', action='store_true',
               default='print available filter specs and filter help and exit')

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


if __name__ == '__main__':
    p = build_parser()
    args = p.parse_args()
    if args.start_after != STATES[0]:
        if args.output_dir is None:
            print('If the --start_after option is used, an existing output directory must be provided (-od).')
            exit()
            args.__dict__['fastq_pref'] = args.output_dir  # ignoring input folder
    else:
        if args.fastq_prefix is None:
            print('If the --start_after option is not used, an input fastq prefix/folder (--fastq_prefix).')
            exit()

    p, s = os.path.split(args.fastq_prefix)
    args.__dict__['fastq_path'] = p
    args.__dict__['fastq_pref'] = s
    args.__dict__['fastq_prefix'] = canonic_path(args.fastq_prefix)
    p, s = os.path.split(args.fastq_prefix)
    if args.barcode_file is None:
        args.__dict__['barcode_file'] = os.sep.join([args.fastq_path,'barcodes'])

    if args.debug is not None:
        nlines, nsamples = args.debug.split(',')
        args.__dict__['db_nlines'] = int(nlines)
        args.__dict__['db_nsamples'] = int(nsamples)

    if args.output_dir is not None:
        args.__dict__['output_dir'] = canonic_path(args.__dict__['output_dir'])
    mh = MainHandler(args, ' '.join(sys.argv))
    mh.execute()