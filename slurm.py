import dill
import sys
import subprocess as sp
from config import *
import os
import shlex as sh
import traceback


def pkl_args(id):
    return '.%s.pkl.args' % id


def pkl_output(id):
    return '.%s.pkl.out' % id


def execute(f, args, kwargs, slurm_spec, interval=.1):
    # prepare python input
    execblob = (f, args, kwargs)
    id = hash(str(execblob))
    with open(pkl_args(id), 'wb') as OUT: dill.dump(execblob, OUT)

    # sbatch script
    tmpscript = '.%s.slurmscript' % id
    with open(tmpscript, 'w') as script:
        script.write('#! /bin/bash \n')
        script.write('%s %s %s\n' % (INTERPRETER, __file__, id))

    # execute
    opts = ' '.join('--%s=%s' % (k,str(v)) for k,v in slurm_spec.items())
    p = sp.Popen(sh.split('sbatch %s %s' % (opts, tmpscript)), stderr=sp.PIPE, stdout=sp.PIPE)
    out, err = p.communicate()

    # monitor until done
    jid, done = out.decode('utf8').strip().split(' ')[-1], False
    while not done:
        q = sp.Popen("squeue",stdout=sp.PIPE).communicate()[0].decode('utf8').split('\n')
        done = True
        for line in q:
            if line.strip().split(' ')[0] == jid:
                done = False
                break

    # handle function output and cleanup
    with open(pkl_output(id), 'rb') as IN: out, err = dill.load(IN)
    os.remove(tmpscript)
    os.remove(pkl_output(id))
    os.remove(pkl_args(id))
    os.remove('slurm-%s.out' % jid)
    if err is not None: raise Exception(err)
    return out


if __name__ == '__main__':
    id = sys.argv[1]
    try:
        err = None
        with open(pkl_args(id), 'rb') as IN:
            f, args, kwargs = dill.load(IN)
            out = f(*args,**kwargs)
    except Exception:
        out, err = None, traceback.format_exc(10)
    with open(pkl_output(id), 'wb') as OUT: dill.dump((out, err), OUT)