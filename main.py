from config import *
import sys

if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    import subprocess as sp
    import os
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()

import argparse
import csv
import logging as lg
import threading as th
import multiprocessing as mp
import threading
import pickle
import shutil
from collections import Counter

from exporters import *
from filters import *
from manage import WorkManager
from secure_smtp import ThreadedTlsSMTPHandler
from utils import *


class FeatureCollection(OrderedDict):

    def __init__(self):
        super(FeatureCollection, self).__init__()

    def add_feature(self, feat):
        for f in self.values():
            if f.short_name == feat.short_name: raise ValueError()
            if f.name == feat.name: raise ValueError()
        self[feat.name] = feat
        if feat.short_name is None:
            i = 1
            short_names = set(f.short_name for f in self.values())
            while i <= len(feat.name):
                n = feat.name[:i].lower()
                if n in short_names: continue
                feat.short_name = n
                break
        if feat.short_name is None:
            raise ValueError()


class Feature(object):

    def __init__(self, name, type, short_name=None, units=None):
        self.name = name.lower()
        self.short_name = short_name
        self.strtype = type.lower()
        self.type = str if type=='str' else int if type == 'int' else float
        self.units = units
        self.vals = set([])

    def __str__(self):
        return '%s(%s)[%s]:%s' % (self.name, self.short_name, self.units, self.strtype)

    def __repr__(self): return str(self)


class Sample(object):

    def __init__(self, context):
        self.fvals = OrderedDict()
        self.barcode = None
        self.context = context
        self.files = {}

    def file_map(self):
        fs = {
            'fastq': os.sep.join([self.context.fastq_dir, self.base_name() + FASTQ_SUFF]),
            'sam_hdr': self.context.bam_dir + os.sep + self.base_name() + SAM_HDR_SUFF,
            'tmp_bam': self.context.tmp_dir + os.sep + self.base_name() + TMP_BAM_SUFF,
            'bam_f': None,
            'unfiltered_bam': self.context.bam_dir + os.sep + self.base_name() + UNFILTERED_SUFF,
            'bam': self.context.bam_dir + os.sep + self.base_name() + BAM_SUFF,
            'cbw': self.context.bw_dir + os.sep + self.base_name() + '.c.bw',
            'wbw': self.context.bw_dir + os.sep + self.base_name() + '.w.bw',
            'tmp_bed': self.context.tmp_dir + os.sep + self.base_name() + TMP_BED_SUFF,
            'tmp_cnt': self.context.tmp_dir + os.sep + self.base_name() + TMP_CNT_SUFF
        }
        if self.context.a.keep_filtered:
            fs['bam_f'] = self.context.filtered_dir + os.sep + self.base_name() + BAM_SUFF
        return fs

    def base_name(self):
        return '_'.join('%s-%s' % (f.short_name, str(v)) for f, v in self.fvals.items())

    def full_name(self):
        return '_'.join('%s-%s' % (f.name, str(v)) for f, v in self.fvals.items())

    def __repr__(self):
        return self.base_name()

    def __hash__(self):
        return hash(tuple(self.fvals.values()))

    def critical(self, msg, err, countq):
        msg += '\n' + err
        self.context.logq.put((lg.CRITICAL, msg))
        countq.put((self.barcode, countq))
        exit()

    def handle(self, in_files, tts_file, countq):
        self.files['in1'] = in_files[0]
        self.files['in2'] = in_files[1]
        self.files.update(self.file_map())
        c = self.context.w_manager.get_channel()
        cp = self.context.w_manager.get_channel() #for parallel tasks

        # fastq
        args = (self.files, self.context.a.umi_length, self.context.bc_len)
        self.context.w_manager.execute(func=Sample.format_fastq, args=args, c=c)
        out, err = c.get()
        if err is not None:
            msg = 'Error while formatting fastq for sample %s' % self.base_name()
            self.critical(msg, err, countq)
        msg = 'Fastq for sample %s is ready: %s' % (self.base_name(), self.files['fastq'])
        self.context.logq.put((lg.INFO, msg))

        # alignment
        args = (self.files, self.context.a.n_threads, self.context.a.align_index_path, self.context.fpipe)
        self.context.w_manager.execute(func=Sample.make_bam, args=args, c=c,
                                       slurm_spec={'cpus-per-task':self.context.a.n_threads,
                                                   'mem':'8G'})
        if self.context.a.count_index_paths is not None:
            args = (self.files, self.context.a.count_index_paths, self.context.a.n_threads)
            self.context.w_manager.execute(func=Sample.make_bam, args=args, c=cp,
                                           slurm_spec={'cpus-per-task': self.context.a.n_threads,
                                                       'mem':'8G'})
            stats, err = cp.get()
            if err is not None:
                msg = ('Error while counting alignment for sample %s' % self.base_name()) + '\n' + err
            else:
                msg = 'Counted alignments for sample %s' % self.base_name()
            self.context.logq.put((lg.INFO, msg))
            self.context.statq.put((self.base_name(), stats))
        stats, err = c.get()
        if err is not None:
            msg = 'Error while aligning sample %s' % self.base_name()
            self.critical(msg, err, countq)
        msg = 'BAM for sample %s is ready: %s' % (self.base_name(), self.files['bam'])
        self.context.logq.put((lg.INFO, msg))
        self.context.statq.put((self.base_name(), stats))

        # tracks and counts
        self.context.w_manager.execute(func=Sample.make_tracks, args=(self.files,), c=c)
        self.context.w_manager.execute(func=Sample.count, args=(tts_file, self.files), c=cp)
        cnt, err = c.get()
        if err is not None:
            msg = ('Error while making tracks for sample %s' % self.base_name()) + '\n' + err
        else:
            msg = 'BigWig tracks ready for sample %s' % self.base_name()
        self.context.logq.put((lg.INFO, msg))

        cnt, err = cp.get()
        if err is not None:
            msg = 'Error while counting tts in sample %s' % self.base_name()
            self.critical(msg, err, countq)
        ttl = sum(int(val) for val in cnt.values())
        msg = 'Counted tts reads in sample %s, total: %s' % (self.base_name(), ttl)
        self.context.logq.put((lg.INFO, msg))
        stats = {'tts_counted': ttl}
        self.context.statq.put((self.base_name(), stats))
        countq.put((self.barcode, cnt))

    @staticmethod
    def format_fastq(files, umi_len, bc_len):
        import shlex as sh
        pre1, pre2 = files['in1'], files['in2']
        if os.path.isfile(pre1) and os.path.isfile(pre2):
            cat = sp.Popen(['cat', pre1, pre2], stdout=sp.PIPE)
        elif os.path.isfile(pre1):
            cat = sp.Popen(['cat', pre1], stdout=sp.PIPE)
        elif os.path.isfile(pre2):
            cat = sp.Popen(['cat', pre2], stdout=sp.PIPE)
        else:
            cat = sp.Popen(['cat'], stdin=open(os.devnull), stdout=sp.PIPE)
        awk = sp.Popen(sh.split('''awk -F "\\t" '{print "@umi:"substr($4,%i,%i)"\\n"$3"\\n+\\n"$7}' '''
                                % (bc_len+1, umi_len)), stdin=cat.stdout, stdout=sp.PIPE)
        gzip = sp.Popen(['gzip'], stdin=awk.stdout, stdout=open(files['fastq'], 'wb'))
        gzip.wait()
        if os.path.isfile(pre1): os.remove(pre1)
        if os.path.isfile(pre2): os.remove(pre2)
        return None

    @staticmethod
    def alignment_count(files, genome_indices, n_threads):
        # align, and parse statistics
        # bowtie2 --local -p 4 -U {fastq.gz} -x {index} 2> {stats} >/dev/null
        import shlex as sh
        from utils import parse_bowtie_stats
        stats = {}
        for genome, index in genome_indices:
            bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (EXEC['BOWTIE'], n_threads, files['fastq'], index)),
                          stderr=sp.PIPE, stdout=open(os.devnull, 'w'))
            tmp_stats = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
            for k, v in tmp_stats:
                stats[genome+'-'+k] = v
        return stats

    @staticmethod
    def make_bam(files, n_threads, genome_index, fpipe):
        import shlex as sh
        from utils import parse_bowtie_stats
        bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (EXEC['BOWTIE'], n_threads, files['fastq'], genome_index)),
                      stdout=sp.PIPE, stderr=sp.PIPE)
        awkcmd = ''.join(("""awk '{if (substr($1,1,1) == "@" && substr($2,1,2) == "SN")""",
                          """{print $0 > "%s";} print; }' """)) % files['sam_hdr']
        geth = sp.Popen(sh.split(awkcmd), stdin=bt.stdout, stdout=sp.PIPE)
        st = sp.Popen(sh.split('samtools view -b -o %s' % files['tmp_bam']), stdin=geth.stdout)
        st.wait()
        if 'unaligned_bam' in files:
            cmd = 'samtools view -f4 -b %s -o %s' % (files['tmp_bam'], files['unaligned_bam'])
            naligned = sp.Popen(sh.split(cmd), stdout=sp.PIPE)
            naligned.wait()
        aligned = sp.Popen(sh.split('samtools view -F4 -b %s' % files['tmp_bam']), stdout=sp.PIPE)
        sort = sp.Popen(sh.split('samtools sort -o %s -T %s' % (files['unfiltered_bam'], files['tmp_bam'])),
                        stdin=aligned.stdout)
        sort.wait()
        os.remove(files['tmp_bam'])
        stats = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
        n = fpipe.filter(files['unfiltered_bam'], files['bam'], files['sam_hdr'], files['bam_f'])
        stats['passed_filter'] = n
        os.remove(files['unfiltered_bam'])
        return stats

    @staticmethod
    def make_tracks(files):
        import shlex as sh

        def handle_strand(char, fout, negate=False):
            bedcmd = "bedtools genomecov -ibam %s -g %s -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam'], SCER_GENOME_LENGTH_PATH, STRANDS[char])), stdout=sp.PIPE)
            if negate:
                bed = sp.Popen(['awk', '{print $1,$2,$3,"-"$4;}'], stdin=bed.stdout, stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(files['tmp_bed'], 'w'))
            sbed.wait()
            bw = sp.Popen([EXEC['BG2W'], files['tmp_bed'], SCER_GENOME_LENGTH_PATH, fout])
            bw.wait()
            os.remove(files['tmp_bed'])

        handle_strand('w', files['wbw'])
        handle_strand('c', files['cbw'], True)

    @staticmethod
    def count(tts_bed, files):
        import shlex as sh
        from collections import OrderedDict
        cnt = sp.Popen(sh.split('bedtools coverage -s -counts -a %s -b %s' % (tts_bed, files['bam'])),
                       stdout=open(files['tmp_cnt'], 'wb'))
        cnt_dict = OrderedDict()
        with open(tts_bed) as ttsf:
            for line in ttsf:
                cnt_dict[line.split('\t')[3]] = 0
        cnt.wait()
        with open(files['tmp_cnt']) as tmp:
            for line in tmp:
                if not line: continue
                sline = line.split('\t')
                cnt_dict[sline[3].strip()] = sline[6].strip()
        os.remove(files['tmp_cnt'])
        return cnt_dict


class ExperimentHandler(object):

    def __init__(self, argobj, cmdline):
        self.a = argobj

        self.logfile, self.logq,  self.logger = self.setup_log()
        self.log(lg.INFO, 'commandline: %s' % cmdline)

        if self.a.debug:
            self.log(lg.INFO, '=================== DEBUG MODE (%s) ===================' % self.a.debug)
        te = check_third_party()
        for t, e in te.items():
            if e is not None:
                self.log(lg.ERROR, '%s: %s' % (t,e))
            else:
                self.log(lg.DEBUG, '%s: %s' % (t, 'OK'))

        self.w_manager = WorkManager(max_w=self.a.max_workers,
                                     default_slurm_spec={'cpus-per-task':2, 'mem':'4G'})

        self.bc_len, self.samples, self.features = self.parse_sample_db()

        self.setup_output()

        self.statq, self.stat_thread = self.setup_stats()

        self.tts_bed_path, self.tts_accs = self.build_tts_file()

        sfname = self.a.output_dir + os.sep + 'sample_db.csv'
        if not os.path.isfile(sfname): shutil.copy(self.a.sample_db, sfname)

        self.fpipe = build_filter_schemes('filter:'+self.a.filter)['filter']
        self.log(lg.INFO, 'Filters:\n' + str(self.fpipe))
        self.exporters = exporters_from_string(self.a.exporters, self.a.output_dir)

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
        if self.a.debug is not None: ch.setLevel(lg.DEBUG)
        if self.a.user_emails is not None:
            mailh = ThreadedTlsSMTPHandler(mailhost=('smtp.gmail.com', 587),
                                           fromaddr='transeq.pipeline@google.com',
                                           toaddrs=self.a.user_emails.split(','),
                                           credentials=('transeq.pipeline', 'transeq1234'),
                                           subject='TRANSEQ pipeline message')
            mailh.setLevel(lg.CRITICAL)
            logger.addHandler(mailh)
        logger.addHandler(fh)
        logger.addHandler(ch)

        def log(q, logger):
            for lvl, msg in iter(q.get, None): logger.log(lvl, msg)

        logq = mp.Queue()
        lt = th.Thread(target=log, args=(logq, logger))
        lt.daemon = True
        lt.start()
        return logfile, logq, lt

    def log(self, lvl, msg):
        self.logq.put((lvl,msg))

    def update_stats(self, sq):
        stats = OrderedDict()
        stat_order = []
        for s in self.samples.values(): stats[s.base_name()] = Counter()
        stats[NO_BC_NAME] = Counter()
        for sname, counter in iter(sq.get, None):
            for k in counter.keys():
                if k not in stat_order: stat_order.append(k)
            stats[sname] += counter
            self.write_stats_file(stats, stat_order)
        sq.put((stats, stat_order))

    def setup_stats(self):
        sq = mp.Queue()
        st = th.Thread(target=self.update_stats, args=(sq,))
        st.daemon = True
        st.start()
        return sq, st

    def build_tts_file(self):
        tts_bed_name = self.tmp_dir + os.sep + 'tts.tmp.bed'
        tts_bed = open(tts_bed_name, 'w')
        chrlens = chr_lengths()
        tts_accs = []
        for line in open(self.a.tts_file):
            acc, chr, orf_start, orf_end, tts = line.strip().split('\t')
            if chr == 'chrXVII': chr = 'chrM'
            orf_start, orf_end = int(orf_start), int(orf_end)
            strand = orf_start < orf_end
            tts_annot_problem = tts == 'NaN' or \
                                (strand and int(tts) < orf_end) or \
                                (not strand and int(tts) > orf_end)
            tts = orf_end if tts_annot_problem else int(tts)
            w_fr = max(0, tts + ((-1) ** (1 - strand)) * self.a.count_window[1 - strand])
            w_to = min(chrlens[chr], tts + ((-1) ** (1 - strand)) * self.a.count_window[strand])
            if not self.a.dont_bound_start:
                if strand:
                    w_fr = max(w_fr, orf_start)
                else:
                    w_to = min(w_to, orf_start)
            tts_bed.write('\t'.join([chr, str(w_fr), str(w_to), acc, '1', '+' if strand else '-']) + '\n')
            tts_accs.append(acc)
        tts_bed.close()
        return tts_bed_name, tts_accs

    def write_stats_file(self, stats, stat_order):
        with open(self.a.output_dir+os.path.sep+'stats.csv','w') as fout:
            wrtr = csv.DictWriter(f=fout, fieldnames=['sample']+stat_order)
            wrtr.writeheader()
            for s, cntr in stats.items():
                row = dict(cntr)
                row.update({'sample':s})
                wrtr.writerow(row)

    def execute(self):
        self.log(lg.INFO, 'Re-compiling fastq files...')
        self.collect_input_fastqs()
        bcout = None
        if self.a.keep_nobarcode:
            bcout = self.fastq_dir + os.sep + NO_BC_NAME + '.fastq.gz'
        self.split_barcodes(no_bc=bcout)

        cq = self.w_manager.get_channel()
        for bc, sample in self.samples.items():
            in_files = [self.tmp_dir + os.sep + sample.base_name() + '-1',
                        self.tmp_dir + os.sep + sample.base_name() + '-2']
            args = (in_files, self.tts_bed_path, cq)
            sample.worker = threading.Thread(target=sample.handle, args=args)
            sample.worker.start()

        tts_counters = {}
        while True:
            bc, counter = cq.get()
            tts_counters[bc] = counter
            if len(tts_counters) == len(self.samples): break

        if self.a.no_hub: self.build_hub()

        self.aftermath(tts_counters)

    def parse_sample_db(self):

        def parse_features(hdr):
            feat_pat = re.compile('\s*(?P<name>\w+)\s*(?:\((?P<short_name>\w+)\))?'
                                  '\s*:(?P<type>\w+)(:?\[(?P<units>\w+)\])?')
            hdr = hdr.split(SAMPLEDB_DELIM)
            features = FeatureCollection()
            f_pos_map = {}
            for i, f in enumerate(hdr):
                if i == 0:
                    assert hdr[0] == 'barcode', 'first column in sample db needs to be the "barcode" column'
                elif f.startswith('#'):
                    msg = ("ignoring column %s in sample db" % f)
                    self.log(lg.INFO, msg)
                else:
                    m = re.match(feat_pat, f)
                    if m is None:
                        msg = ("couldn't understand feature '%s' in sample_db file, format should be: "
                               "<name>(<short_name>):(str|int|float)[units] (short_name and units are optional) or "
                               "column is ignored if it starts with '#'" % f)
                        self.log(lg.CRITICAL, msg)
                        raise (ValueError(msg))
                    try:
                        f_pos_map[i] = Feature(**m.groupdict())
                    except ValueError:
                        snames = '\n'.join(f.short_name for f in features.values)
                        msg = ("features must have distinct names and short_names - %s appears at least twice (or "
                               "its short_name matched a previous generated short name):\n%s" % f, snames)
                        self.log(lg.CRITICAL, msg)
                        raise (ValueError(msg))
                    features.add_feature(f_pos_map[i])
            return features, f_pos_map

        def parse_samples(file, f_pos_map):
            b2s, bc_len = OrderedDict(), None
            for i, line in enumerate(file):
                if line.strip()[0] == '#': continue  # comment
                if self.a.debug is not None:
                    if i >= self.a.db_nsamples: break  # limit number of samples
                sample = Sample(context=self)
                for j, val in enumerate(line.strip().split(SAMPLEDB_DELIM)):
                    val = val.strip()
                    if j == 0:
                        if i == 0:
                            bc_len = len(val)  # first barcode
                        elif bc_len != len(val):
                            msg = "barcode %s has a different length" % val
                            self.log(lg.CRITICAL, msg)
                            raise (TypeError(msg))
                        if val in b2s:
                            msg = "barcode %s is not unique" % val
                            self.log(lg.CRITICAL, msg)
                            raise (TypeError(msg))
                        sample.barcode = val
                    elif j in f_pos_map:
                        f = f_pos_map[j]
                        try:
                            v = f.type(val)
                        except ValueError:
                            msg = ("couldn't cast value %s in sample %i, feature '%s' to "
                                   "given type - %s." % (val, i + 1, f.name, f.strtype))
                            self.log(lg.CRITICAL, msg)
                            raise (ValueError(msg))
                        f.vals.add(v)
                        sample.fvals[f] = v
                    if hash(sample) in [hash(s) for s in b2s.values()]:
                        msg = "2 samples (or more) seem to be identical - %s" % sample
                        self.log(lg.CRITICAL, msg)
                        raise (TypeError(msg))
                b2s[sample.barcode] = sample
            return b2s, bc_len

        sdb = open(self.a.sample_db)
        exp = re.match('.*experiment.*:\s+(\w+)', sdb.readline())
        if exp is None:
            msg = 'barcodes file should contain a header with experiment name: ' \
                  'experiment: <expname>'
            self.log(lg.CRITICAL, msg)
            raise ValueError(msg)
        self.user = getpass.getuser()
        self.exp = exp.group(1)
        msg = 'user: %s, experiment: %s' % (self.user, self.exp)
        self.log(lg.INFO, msg)
        features, f_pos_map = parse_features(sdb.readline())
        b2s, bc_len = parse_samples(sdb, f_pos_map)
        sdb.close()
        msg = '\n'.join(['barcodes:'] + ['%s -> %s' % (b, s.base_name()) for b, s in b2s.items()])
        self.log(lg.DEBUG, msg)
        self.log(lg.INFO, 'found %i samples.' % len(b2s))
        msg = 'features:\n' + '\n'.join('%s: %s' % (str(f), ','.join(str(x) for x in f.vals)) for f in features.values())
        self.log(lg.DEBUG, msg)

        return bc_len, b2s, features

    def collect_input_fastqs(self):
        files = {}
        for fn in os.listdir(self.a.fastq_path):
            if not fn.endswith('fastq.gz'): continue
            if not fn.startswith(self.a.fastq_pref): continue
            parts = re.split('_R\d',fn)
            if len(parts) == 2:
                path = self.a.fastq_path + os.sep + fn
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
        self.log(lg.INFO, msg)
        if not files:
            msg = "could not find R1/R2 fastq.gz pairs in given folder: %s" % self.a.fastq_path
            self.log(lg.CRITICAL, msg)
            raise IOError(msg)
        self.input_files = files

    def setup_output(self):
        if self.a.output_dir is None:
            folder = canonic_path(DATA_PATH) + os.sep + self.user
            self.dir_and_log(folder)
            folder += os.sep + self.exp
            self.dir_and_log(folder)
            folder += os.sep + datetime.datetime.now().strftime("%d-%m-%y")
            self.dir_and_log(folder)
            if self.a.debug is not None:
                fname = 'debug'
            else:
                i = [int(x) for x in os.listdir(folder) if isint(x)]
                if not i: i = 1
                else: i = max(i) + 1
                fname = str(i)
            folder += os.sep + fname
            self.a.output_dir = folder

        d = self.a.output_dir
        self.tmp_dir = d + os.sep + TMP_DIR
        self.fastq_dir = d + os.sep + FASTQ_DIR
        self.bw_dir = d + os.sep + BW_DIR
        self.bam_dir = d + os.sep + BAM_DIR

        if os.path.isdir(self.tmp_dir): shutil.rmtree(self.tmp_dir)
        if os.path.islink(self.bw_dir): os.remove(self.bw_dir)
        elif os.path.isdir(self.bw_dir): shutil.rmtree(self.bw_dir)

        # assuming all folder structure exists if check passes
        self.dir_and_log(d, lg.INFO)
        self.dir_and_log(self.fastq_dir)
        self.dir_and_log(self.bam_dir)
        self.dir_and_log(self.bw_dir)
        self.dir_and_log(self.tmp_dir)
        if self.a.keep_filtered:
            self.filtered_dir = d + os.sep + FILTERED_DIR
            self.dir_and_log(self.filtered_dir)
        if self.a.keep_unaligned:
            self.unaligned_dir = d + os.sep + UNALIGNED_DIR
            self.dir_and_log(self.unaligned_dir)

        if self.a.no_hub:
            if self.a.hub_name is None: self.a.hub_name = self.exp
            self.www_rel = os.sep.join([getpass.getuser(), self.a.hub_name])
            d = os.sep.join([WWW_PATH, getpass.getuser()])
            self.dir_and_log(d, lg.INFO,chto='777')
            d = os.sep.join([d, self.a.hub_name])
            if os.path.isdir(d):
                self.log(lg.DEBUG, 'Removing old folder')
                shutil.rmtree(d)
            self.dir_and_log(d, lg.INFO, chto='777')
            self.www_path = d

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
            cntr = Counter()
            with open(bc1) as IN:
                for line in IN:
                    sample, cnt = line.strip().split(' ')
                    sample = sample[:-2]
                    if sample == NO_BC_NAME: continue  # only from nobc_counts-2
                    c = Counter()
                    c[stat] = int(cnt)
                    cntr[sample] += int(cnt)
                    self.statq.put((sample, c))
            os.remove(bc1)
            with open(bc2) as IN:
                for line in IN:
                    sample, cnt = line.strip().split(' ')
                    c = Counter()
                    c[stat] = int(cnt)
                    cntr[sample[:-2]] += int(cnt)
                    self.statq.put((sample[:-2], c))
            os.remove(bc2)
            for s in self.samples.values():
                if s.base_name() not in cntr:
                    cntr[s.base_name()] += 0
            msg = '\n'.join(['%s: %i' % (s, c) for s, c in cntr.items()])
            self.log(lg.CRITICAL, 'read counts:\n' + msg)

        hb = {}
        for b,s in self.samples.items():
            hb.update({eb:s.base_name() for eb in hamming_ball(b, self.a.hamming_distance)})

        awk1p, cnt1 = compile_awk("1", {b: s.base_name() for b,s in self.samples.items()})
        awk2p, cnt2 = compile_awk("2", hb)
        outf = open(os.devnull, 'w') if no_bc is None else open(no_bc, 'wb')
        for r1, r2 in self.input_files:
            msg = 'splitting files:\n%s\n%s' % (os.path.split(r1)[1],os.path.split(r2)[1])
            self.log(lg.INFO, msg)
            paste1 = sp.Popen('paste <(zcat %s) <(zcat %s)' % (r1,r2), stdout=sp.PIPE,
                              shell=True, executable='/bin/bash')
            awkin = sp.Popen(sh.split('paste - - - -'), stdin=paste1.stdout, stdout=sp.PIPE)
            if self.a.debug: # only a small subset of reads
                nlines = round(self.a.db_nlines/4)*4  # making sure it's in fastq units
                awkin = sp.Popen(sh.split('head -%i' % nlines), stdin=awkin.stdout, stdout=sp.PIPE)
            awk1 = sp.Popen(sh.split('awk -F "\\t" ' + awk1p), stdin=awkin.stdout, stdout=sp.PIPE)
            awk2 = sp.Popen(sh.split('awk -F "\\t" ' + awk2p), stdin=awk1.stdout, stdout=sp.PIPE)
            awkcmd = """awk -F "\\t" '{print $1"\\n"$3"\\n"$5"\\n"$7; print $2"\\n"4"\\n"$6"\\n"$8;}' """
            wfastq = sp.Popen(sh.split(awkcmd), stdin=awk2.stdout, stdout=sp.PIPE)
            gzip = sp.Popen(['gzip'], stdin=wfastq.stdout, stdout=outf)
            wfastq.wait() # to prevent data interleaving
        gzip.wait()
        self.log(lg.INFO, 'Barcode splitting finished.')

        merge_statistics(cnt1, cnt2)

    def build_hub(self):
        self.log(lg.INFO, 'Generating hub...')
        sacpath = self.www_path + os.path.sep + 'sacCer3'
        if not os.path.exists(sacpath): os.mkdir(sacpath)
        hubfile = open(self.www_path + os.path.sep + 'hub.txt', 'w')
        hubfile.write('\n'.join(["hub %s" % self.a.hub_name,
                                 "shortLabel %s" % self.a.hub_name,
                                 "longLabel %s(%s)" % (self.a.hub_name, self.exp),
                                 "genomesFile genomes.txt",
                                 "email %s" % self.a.hub_email]))
        genomesfile = open(self.www_path + os.path.sep + 'genomes.txt', 'w')
        genomesfile.write("genome sacCer3\n"
                          "trackDb sacCer3/trackDB.txt")
        trackfile = open(sacpath + os.path.sep + 'trackDB.txt', 'w')
        for s in self.samples.values():
            wurl = os.sep.join([URL_BASE, self.www_rel, BW_DIR, s.base_name()+'.w.bw'])
            curl = os.sep.join([URL_BASE, self.www_rel, BW_DIR, s.base_name()+'.c.bw'])
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
        mainurl = os.sep.join([URL_BASE, self.www_rel, 'hub.txt'])
        new_path = self.www_path + os.sep + BW_DIR
        shutil.move(self.bw_dir, self.www_path)
        os.symlink(new_path, self.bw_dir, target_is_directory=True)
        sp.call('chmod -R 777 %s' % self.www_path, shell=True)
        msg = 'Transferred bigwig files to %s (link available in %s as well)' % (new_path, self.a.output_dir)
        self.log(lg.DEBUG, msg)
        self.log(lg.CRITICAL, 'Hub available at %s' % mainurl)

    def export(self, tts_cnts, out_stats, stat_order):
        stats = OrderedDict()
        cnts = OrderedDict()
        for s in self.samples.values():
            stats[s] = out_stats[s.base_name()]
            cnts[s] = tts_cnts[s.barcode]
        all = [('stats', stats, stat_order),
               ('tts', cnts, self.tts_accs)]

        for e in self.exporters:
            fs = e.export(self.features.values(), self.samples.values(), all)
            for f in fs:
                self.log(lg.INFO, 'Exported data to file: %s' % (self.a.output_dir+os.sep+f,))
                if self.a.export_path is not None:
                    target = self.a.export_path + os.sep + self.exp + '-' + f
                    shutil.copy(self.a.output_dir + os.sep + f, target)
                    self.log(lg.DEBUG, 'Copied data to: %s' % target)

    def aftermath(self, tts_counters):
        # remove temp folder
        # modify file permissions for the entire tree
        # make everything read only (optional?)
        # store pipeline code?
        # merge and report statistics
        # merge data to single (usable) files
        if self.a.debug is None:
            shutil.rmtree(self.tmp_dir)

        pickle.dump(self.a, open(self.a.output_dir + os.sep + 'args.pkl', 'wb'))

        self.log(lg.CRITICAL, 'All done.')
        shutil.copy(self.logfile, self.a.output_dir + os.sep + 'full.log')
        self.logq.put(None)
        self.logger.join()
        self.statq.put(None)
        all_stats, sord = self.statq.get()
        self.export(tts_counters, all_stats, sord)
        self.stat_thread.join()
        self.w_manager.close()

    def dir_and_log(self, path, level=lg.DEBUG, chto='770'):
        self.log(level, 'creating folder %s' % path)
        create_dir(path, chto)


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('fastq_prefix', type=str, default=None,
                   help='path to a prefix of fastq files (R1 & R2) containing the transeq data.'
                        'This can be a folder (must end with "/"), in which case all R1/R2 pairs'
                        'in the folder are considered, or a "path/to/files/prefix", in which case '
                        'all files in the path with the prefix are considered')
    g.add_argument('--max_workers', '-mw', type=int, default=100,
                   help='maximal number of simultaneous working processes in this pipeline')
    g = p.add_argument_group('Output')
    g.add_argument('--output_dir', '-od', default=None, type=str,
                   help='path to the folder in which most files are written. '
                        'If not given, the date and info from barcode file are used to '
                        'generate a new folder in %s' % DATA_PATH)
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
    g.add_argument('--align_index_path', '-aip', type=str, default='/cs/wetlab/genomics/scer/bowtie/sacCer3',
                   help='path prefix to bowtie genome index for alignment, default is'
                        '/cs/wetlab/genomics/scer/bowtie/sacCer3 ')
    g.add_argument('--count_index_paths', '-cip', type=str, default=None,
                   help='comma separated list of <genome_name>:<index_path> for counting alignment events. '
                        'For example: '
                        '"klac:/cs/wetlab/genomics/klac/bowtie/genome,human:/cs/wetlab/genomics/human/bowtie/genome"')
    g.add_argument('--n_threads', '-an', type=int, default=4,
                   help='number of threads used for alignment per bowtie instance')
    g.add_argument('--keep_unaligned', '-ku', action='store_true',
                   help='if set, unaligned reads are written to '
                        'output_folder/%s/<sample_name>.bam' % UNALIGNED_DIR)

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
                        'output_folder/%s/<sample_name>.bam' % FILTERED_DIR)
    g.add_argument('--filter_specs', '-fh', action='store_true',
                   help='print available filters, filter help and exit')

    g = p.add_argument_group('Tracks')
    g.add_argument('--no_hub', '-mh', action='store_false',
                   help='prevent the pipeline from genrating a browser hub in the default www folder')
    g.add_argument('--hub_email', '-he', action='store', default='noemail@nodomain.com',
                   help='the contact email for the generated hub')
    g.add_argument('--hub_name', '-hn', action='store', default=None,
                   help='the directory in which the hub is generated. If not given, '
                        'experiment name is used.')

    g = p.add_argument_group('Count')
    g.add_argument('--tts_file', '-tf', default=TTS_MAP,
                   help='annotations for counting. Expected format is a tab delimited file with "chr", "ACC", "start",'
                        '"end", and "TTS" columns. default is found at %s' % TTS_MAP)
    g.add_argument('--keep_tts_bed', '-ktb', action='store_true',
                   help='whether the TTS window definition bed file should be kept or dicarded')
    g.add_argument('--dont_bound_start', '-dbs', action='store_false',
                   help='do not limit counting window orf start (i.e. if longer than ORF, trim counting window)')
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

    if args.count_index_paths is not None:
        args.__dict__['count_index_paths'] = [pair.split(':') for pair in args.count_index_paths.split(',')]

    args.__dict__['count_window'] = [int(x) for x in args.count_window[1:-1].split(',')]

    return args


if __name__ == '__main__':
    p = build_parser()
    a = parse_args(p)
    mh = ExperimentHandler(a, ' '.join(sys.argv)) #also print git current version from ./.git/log/HEAD,last row, 2nd column
    mh.execute()