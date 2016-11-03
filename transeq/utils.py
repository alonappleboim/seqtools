import datetime
import getpass
import os
import re
import shlex as sh
import subprocess as sp

from common.config import *


def chr_lengths():
    cl = {}
    with open(SCER_GENOME_LENGTH_PATH) as cl_file:
        for line in cl_file:
            c, l = line.strip().split('\t')
            cl[c] = int(l)
    return cl


def buff_lines(buffer, buf_size=1024):
    prev = ''
    while True:
        data = buffer.read(buf_size)
        if not data: break
        else:
            lines = ''.join(data.decode('utf-8')).split('\n')
            lines[0] = prev + lines[0]
            prev = lines.pop(-1)
            for line in lines: yield line
    yield prev


def create_dir(path, chto='770'):
    if not os.path.exists(path):
        os.mkdir(path)
        sp.Popen(sh.split('chmod -R %s %s' % (chto, path)), stderr=sp.PIPE).communicate()


def get_logfile():
    now = datetime.datetime.now()
    p = LOG_PATH + os.sep + str(now.year)
    create_dir(p)
    sp.Popen(sh.split('chmod -R 770 %s' % p), stderr=sp.PIPE).communicate()
    p += os.sep + str(now.month)
    create_dir(p)
    sp.Popen(sh.split('chmod -R 770 %s' % p), stderr=sp.PIPE).communicate()
    fname = '%s-%i-%ih%im%s' % (getpass.getuser(), now.day, now.hour, now.minute, now.second)
    return p + os.sep + fname


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


def parse_bowtie_stats(bt_stat_output):
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
    BT_STAT_MAP = {
        2: 'total',
        3: 'no-align',
        4: 'unique-align',
        5: 'multiple-align',
    }
    nlines, stats = 0, {}
    for line in bt_stat_output:
        m = re.match('\s*(\d+).*', line)
        if m is not None:
            if int(m.group(1)) == 0:
                return {s:0 for s in BT_STAT_MAP.values()} # no reads...
            nlines += 1
            if nlines in BT_STAT_MAP.keys():
                stats[BT_STAT_MAP[nlines]] = int(m.group(1))
    return stats


def canonic_path(fname):
    return os.path.abspath(os.path.expanduser(fname))


if __name__ == '__main__':
    for i,line in enumerate(buff_lines(open('README','rb'), 256)): print(line)
    print(i)


