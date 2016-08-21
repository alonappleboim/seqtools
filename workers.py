import time
import os
import subprocess as sp
import multiprocessing as mp
import shlex as sh
import re
import threading
import traceback

from config import *

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
        for _ in range(RETRIALS):
            try:
                e = None
                info = f(**kwargs)
                break  # ok!
            except Exception:
                time.sleep(RETRY_INTERVAL)
                e = traceback.format_exc(10)
        if e is not None: comq.put((token, e, None))
        else: comq.put((token, None, info))


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


def align(files, bowtie_exec, n_threads, scer, klac):

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

    def klac_count():
        # align, and parse statistics
        # bowtie2 --local -p 4 -U {fastq.gz} -x {index} 2> {stats} >/dev/null
        klacstats = files['klac_align_stats']
        bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (bowtie_exec, n_threads, files['fastq'], klac)),
                      stderr=sp.PIPE, stdout=open(os.devnull, 'w'))
        Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
        with open(klacstats, 'w') as S:
            S.write('\n'.join(['total\t%s' % Ns[0],
                               'no-align\t%s' % Ns[1],
                               'unique-align\t%s' % Ns[2],
                               'multiple-align\t%s' % Ns[3]]))

    if klac is not None:
        ct = threading.Thread(target=klac_count)
        ct.start()
        bt = sp.Popen(sh.split('%s --local -p %i -U %s -x %s' % (bowtie_exec, n_threads, files['fastq'], scer)),
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
    sort = sp.Popen(sh.split('samtools sort -o %s -T %s' % (files['bam'], files['tmp_bam'])),
                    stdin=aligned.stdout)
    sort.wait()
    ind = sp.Popen(sh.split('samtools index %s' % files['bam']))
    ind.wait()
    os.remove(files['tmp_bam'])
    Ns = parse_bowtie_stats(''.join(bt.stderr.read().decode('utf8')).split('\n'))
    with open(files['align_stats'], 'w') as S:
        S.write('\n'.join(['total\t%s' % Ns[0],
                           'no-align\t%s' % Ns[1],
                           'unique-align\t%s' % Ns[2],
                           'multiple-align\t%s' % Ns[3]]))
    if klac is not None:
        ct.join()  # making sure lactis alignment is also complete
    return None


def filter_bam(files, fscheme):
    return fscheme.filter(files['bam_in'], files['bam_out'], files['sam_hdr'])


def analyze_bam(files, analyzer):
    return analyzer.analyze(files)
