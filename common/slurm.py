'''
This module abstracts usage of a slurm cluster. Usage:
---
import slurm
def f(x,y): sum(i in range(x*y))
print(55==slurm.execute(f, args=(2,5))) # True!
---
It's also possible to pass key-word arguments to f.
CRITICAL: if f uses any non built-in modules, they must be imported within f
''' #TODO: change this...

import os
import sys

project_dir = os.path.sep.join(sys.modules[__name__].__file__.split(os.path.sep)[:-2])
sys.path.append(project_dir)
from common.config import *

import time
import shlex as sh
import subprocess as sp
import sys
import traceback

import dill


RESPONSE_TIME = .25 # sec

def pkl_args(id):
    return '.%s.pkl.args' % id


def pkl_output(id):
    return '.%s.pkl.out' % id


def execute(f, args, kwargs, slurm_spec, interval=.1):
    # prepare python input
    execblob = (f, args, kwargs)
    id = hash(str(execblob))
    id = id ** 2
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

    # monitor until done and file is available
    jid, done = out.decode('utf8').strip().split(' ')[-1], False
    while not done:
        time.sleep(RESPONSE_TIME)
        q = sp.Popen("squeue",stdout=sp.PIPE).communicate()[0].decode('utf8').split('\n')
        done = True
        for line in q:
            if line.strip().split(' ')[0] == jid:
                done = False
                break
    while not os.path.isfile(pkl_output(id)): time.sleep(RESPONSE_TIME)

    # handle function output and cleanup
    try:
        with open(pkl_output(id), 'rb') as IN: out, err = dill.load(IN)
    except IOError:
        with open('slurm-%s.out' % jid) as LOG: err = 'slurm error: \n %s' % LOG.read()
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
