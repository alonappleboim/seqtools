import os
import abc
from collections import OrderedDict
from common.config import *
import logging as lg
import subprocess as sp

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
    __class__ = abc.ABCMeta

    def __init__(self, context):
        self.fvals = OrderedDict()
        self.barcode = None
        self.context = context
        self.files = {}

    @abc.abstractmethod
    def file_map(self):
         pass

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
        countq.put((self.barcode, {}))
        exit()

    @abc.abstractmethod
    def handle(self, in_files, *args, **kwargs):
        pass

class TranseqSample(Sample):

    def __init__(self, context):
        super(TranseqSample, self).__init__(context)

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
        from transeq.utils import parse_bowtie_stats
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
        from transeq.utils import parse_bowtie_stats
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
            bed = sp.Popen(sh.split(bedcmd % (files['bam'], COMMON_GENOMES['SCER']['chrlens'], STRANDS[char])), stdout=sp.PIPE)
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
