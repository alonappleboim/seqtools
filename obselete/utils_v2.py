import datetime
import getpass
import os
import re
import shlex as sh
import subprocess as sp

from common.config import *

ERROR = -1
BEGIN = 0
FASTQ = 1
ALIGN = 2
COUNT = 3

def create_dir(path):
    if not os.path.exists(path): os.mkdir(path)


def buffered_lines(buffer, buf_size=1024):
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


def parse_bowtie_stats(bytestream):
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
    for i, line in enumerate(''.join(bytestream.read().decode('utf8')).split('\n')):
        m = re.match('\s*(\d+).*', line)
        if m is not None:
            if int(m.group(1)) == 0: return [0, 0, 0, 0]  # no reads...
            nlines += 1
            if nlines in [2, 3, 4, 5]:
                Ns.append(m.group(1))
    return {'no-align':Ns[1], 'unique-align':Ns[2], 'multiple-align':Ns[3]}


def get_logfile(string=''):
    now = datetime.datetime.now()
    p = LOG_PATH + os.sep + str(now.year)
    create_dir(p)
    p += os.sep + str(now.month)
    create_dir(p)
    fname = '%s-%i-%ih%im%s-%s' % (getpass.getuser(), now.day, now.hour, now.minute, now.second, string)
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


if __name__ == '__main__':
    for i,line in enumerate(buff_lines(open('README','rb'), 256)): print(line)
    print(i)


