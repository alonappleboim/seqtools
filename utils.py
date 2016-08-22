import getpass
import datetime
import os
from config import *


def create_dir(path):
    if not os.path.exists(path): os.mkdir(path)


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


