from collections import OrderedDict
import os
import subprocess as sp
import multiprocessing as mp
import shlex as sh
import re
import threading
import traceback

from config import *
from config import SCER_GENOME_LENGTH_PATH as SGLP

MAX_TRIALS = 3


class WorkManager(object):

    def __init__(self, comq, is_slurm=False):
        self.comq = comq
        self.is_slurm = is_slurm
        self.token = 0

    def get_token(self):
        t = self.token
        self.token += 1
        return t

    def run(self, func, kwargs):
        tok = self.get_token()
        if not self.is_slurm:
            w = mp.Process(target=WorkManager.exec_wrapper,
                           args=(func, kwargs, self.comq, tok))
            w.daemon = True
            w.start()
        else:
            pass

        return tok

    @staticmethod
    def exec_wrapper(f, kwargs, comq, token):
        try:
            info = f(**kwargs)
            comq.put((token, None, info))
        except Exception:
            comq.put((token, traceback.format_exc(10), None))


def format_fastq(files, bc_len, umi_len):
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
                            % (bc_len + 1, umi_len)), stdin=cat.stdout, stdout=sp.PIPE)
    gzip = sp.Popen(['gzip'], stdin=awk.stdout, stdout=open(files['fastq'], 'wb'))
    gzip.wait()
    if os.path.isfile(pre1): os.remove(pre1)
    if os.path.isfile(pre2): os.remove(pre2)
    return None


def make_bam(files, bowtie_exec, n_threads, genome, other, fpipe):

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
                if int(m.group(1)) == 0: return [0, 0, 0, 0]  # no reads...
                nlines += 1
                if nlines in [2, 3, 4, 5]:
                    Ns.append(m.group(1))
        return Ns

    def other_count():
        # align, and parse statistics
        # bowtie2 --local -p 4 -U {fastq.gz} -x {index} 2> {stats} >/dev/null
        stats = files['other_stats']
        bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (bowtie_exec, n_threads, files['fastq'], other)),
                      stderr=sp.PIPE, stdout=open(os.devnull, 'w'))
        Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
        with open(stats, 'w') as S:
            S.write('unique-align\t%s' % Ns[2]) #not thread safe!

    if other is not None:
        ct = threading.Thread(target=other_count)
        ct.start()
    bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (bowtie_exec, n_threads, files['fastq'], genome)),
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
    Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
    with open(files['align_stats'], 'w') as S:
        S.write('\n'.join([#'total\t%s' % Ns[0], no need, already added as n_reads
                           'no-align\t%s' % Ns[1],
                           'unique-align\t%s' % Ns[2],
                           'multiple-align\t%s' % Ns[3]]))
    if other is not None:
        ct.join()  # making sure other alignment is also complete
    n = fpipe.filter(files['unfiltered_bam'], files['bam'], files['sam_hdr'], files['bam_f'])
    os.remove(files['unfiltered_bam'])
    return n


def make_tracks(files):

    def handle_strand(char, fout, negate=False):
        bedcmd = "bedtools genomecov -ibam %s -g %s -bg -strand %s"
        bed = sp.Popen(sh.split(bedcmd % (files['bam'], SGLP, STRANDS[char])), stdout=sp.PIPE)
        if negate:
            bed = sp.Popen(['awk', '{print $1,$2,$3,"-"$4;}'], stdin=bed.stdout, stdout=sp.PIPE)
        sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(files['tmp_bed'], 'w'))
        sbed.wait()
        bw = sp.Popen([BG2W_EXEC, files['tmp_bed'], SGLP, fout])
        bw.wait()
        os.remove(files['tmp_bed'])

    handle_strand('w', files['wbw'])
    handle_strand('c', files['cbw'], True)
    return {}


def count(annot_file, files):
    cnt = sp.Popen(sh.split('bedtools coverage -s -counts -a %s -b %s' % (annot_file, files['bam'])),
                   stdout=open(files['tmp_cnt'],'wb'))
    cnt_dict = OrderedDict()
    with open(annot_file) as ttsf:
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

