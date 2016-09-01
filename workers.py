from collections import OrderedDict
import os
import sys
import subprocess as sp
import multiprocessing as mp
from queue import Empty
import shlex as sh
import re
import threading
import traceback

from config import *
from config import SCER_GENOME_LENGTH_PATH as SGLP

MAX_TRIALS = 3


class WorkManager(object):
    """
    An object that allows one to send parallel tasks without worrying about where and when they are executed.
    To use it:
    wm = WorkManager()
    wm.exec(func1, kwargs=dict())
    d# do stuff without waiting
    wm.exec(func1, kwargs=dict(), report_q=myQ)
    result, err = myQ.get() # waiting for result
    if err is None: handle(result)
    """

    def dispatch(self):
        """
        constantly try and execute new tasks from the work queue until notified by manager, whilst keeping
        number of active workers below set number
        """
        wid, roster, tasks, incoming = 0, {}, [], True
        intercom = self.get_channel('internal')
        while incoming or roster or tasks:
            # executing tasks
            while len(roster) < self.max_w and tasks:
                func, kwargs, repq_name = tasks.pop()
                w = mp.Process(target=self.exec_wrapper,
                               args=(func, kwargs, repq_name, wid, self.use_slurm))
                w.start()
                roster[wid] = w
                wid += 1
            # collect new tasks
            if incoming:
                while True:
                    try:
                        task = self.radio['work'].get(timeout=self.delay)
                        if task is None: incoming = False  # no more incoming tasks
                        else: tasks.append(task)
                    except Empty: break
            # remove completed tasks from roster
            while True:
                try: del roster[intercom.get(timeout=self.delay)]
                except Empty: break

    def __init__(self, use_slurm=False, max_w=sys.maxsize, delay=.1):
        self.manager = mp.Manager()
        self.radio = {}
        self.get_channel('work')
        self.use_slurm = use_slurm
        self.max_w = max_w
        self.delay = delay
        self.dispatcher = threading.Thread(target=self.dispatch)
        self.dispatcher.start()

    def get_channel(self, qname):
        self.radio[qname] = self.manager.Queue()
        return self.radio[qname]

    def close(self):
        self.radio['work'].put(None)

    def execute(self, func, kwargs, repq_name=None):
        """
        execute func with kwargs, and report results to the mp.Queue "report_q"
        """
        self.radio['work'].put((func, kwargs, repq_name))

    def exec_wrapper(self, f, kwargs, repq_name, wid, use_slurm):
        err = None  # benefit of the doubt
        kwargs['radio'] = self.radio
        try:
            if use_slurm: out = WorkManager.slurm_exec(f, kwargs)
            else: out = f(**kwargs)
        except Exception:
            out, err = None, traceback.format_exc(10)
        if repq_name is not None: self.radio[repq_name].put((out, err))
        self.radio['internal'].put(wid)  # notify dispatcher that I'm done

    @staticmethod
    def slurm_exec(): pass


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


def make_bam(files, bowtie_exec, n_threads, scer, klac, fpipe):

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
            S.write('unique-align\t%s' % Ns[2])

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
    if klac is not None:
        ct.join()  # making sure lactis alignment is also complete
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
    cnt = sp.Popen(sh.split('bedtools coverage -counts -a stdin -b %s' % files['bam']),
                   stdin=open(annot_file), stdout=sp.PIPE)

    cnt_dict = OrderedDict()
    with open(annot_file) as ttsf:
        for line in ttsf:
            cnt_dict[line.split('\t')[3]] = 0

    for line in ''.join(cnt.stdout.read().decode('utf-8')).split('\n'):
        if not line: continue
        sline = line.split('\t')
        cnt_dict[sline[3].strip()] = sline[6].strip()
    cnt.wait()

    return cnt_dict


if __name__ == '__main__':
    import time


    def log(lq):
        for msg in iter(lq.get,None):
            print(msg)

    def f(x, radio=None):
        channel = tq.get()
        radio[channel].put(print(x**2))
        time.sleep((x*3)**.5)
        radio[channel].put(print(x**.5))

    wm = WorkManager(False, 2, 1)
    lq = wm.get_channel('log')
    tq = wm.get_channel('tmp')
    lp = mp.Process(target=log, args=(lq,))
    lp.start()
    for i in range(5):
        wm.execute(f, {'x':i})
        tq.put('log')
    wm.close()
