# if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
#     import subprocess as sp
#     import os
#     scriptpath = os.path.abspath(sys.modules[__name__].__file__)
#     sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
#     exit()

import functools
import inspect
import multiprocessing as mp
import sys
import threading
import traceback
from queue import Empty

import dill

from common import slurm


def dict_annotated_function(default_setter=None):
    """
    Decorator to verify parameter types, and optionally set defaults
    default_setter(f) can access and change arguments via: f.arguments['argname'] = argval
    """

    def decorator(func):
        sig = inspect.signature(func)

        types = {}
        for param in sig.parameters.values():
            # Iterate through function's parameters and build the list of
            # arguments types
            if (param.annotation is param.empty or
                    type(param.annotation) is not dict or
                    'type' not in param.annotation or
                    not inspect.isclass(param.annotation['type'])):
                continue
            types[param.name] = param.annotation['type']

        def check_type(sig, arg_name, arg_type, arg_value):
            # Internal function that encapsulates arguments type checking
            if not isinstance(arg_value, arg_type):
                raise ValueError("{func}: wrong type of {arg!r} argument, "
                                 "{exp!r} expected, got {got!r}".
                                 format(func=func.__qualname__, arg=arg_name,
                                        exp=arg_type.__name__, got=type(arg_value).__name__))

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # bind the arguments
            ba = sig.bind(*args, **kwargs)
            if default_setter is not None: default_setter(ba)
            for arg_name, arg in ba.arguments.items():
                # And iterate through the bound arguments
                try: type_ = types[arg_name]
                except KeyError: continue
                else:
                    param = sig.parameters[arg_name]
                    if param.kind == param.VAR_POSITIONAL:
                        # If this parameter is a variable-argument parameter,
                        # then we need to check each of its values
                        for value in arg:
                            check_type(sig, arg_name, type_, value)
                    elif param.kind == param.VAR_KEYWORD:
                        # If this parameter is a variable-keyword-argument parameter:
                        for subname, value in arg.items():
                            check_type(sig, arg_name + ':' + subname, type_, value)
                    else:
                        # And, finally, if this parameter a regular one:
                        check_type(sig, arg_name, type_, arg)

            result = func(*ba.args, **ba.kwargs)

            # The last bit - let's check that the result is correct
            return_type = sig.return_annotation
            if (return_type is not sig.empty and
                    isinstance(return_type, type) and
                    not isinstance(result, return_type)):
                raise ValueError('{func}: wrong return type, {exp} expected, got {got}'. \
                                 format(func=func.__qualname__, exp=return_type.__name__,
                                        got=type(result).__name__))
            return result

        return wrapper

    return decorator


class ExecError(Exception): pass


class WorkManager(object):
    """
    An object that allows one to send parallel tasks without worrying about how, where and when they are executed.
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
        number of active work below set number
        """
        wid, roster, tasks, incoming = 0, {}, [], True
        intercom = self.get_channel()
        while incoming or roster or tasks:
            # executing tasks
            while len(roster) < self.max_w and tasks:
                slurm_spec, func, args, kwargs, c = dill.loads(tasks.pop())
                w = mp.Process(target=self.exec_wrapper,
                               args=(func, args, kwargs, c, intercom, wid, slurm_spec))
                w.start()
                roster[wid] = w
                wid += 1
            # collect new tasks
            if incoming:
                while True:
                    try:
                        task = self.work.get(timeout=self.delay)
                        if task is None: incoming = False  # no more incoming tasks
                        else: tasks.append(task)
                    except Empty: break
            # remove completed tasks from roster
            while True:
                try: del roster[intercom.get(timeout=self.delay)]
                except Empty: break
        self.work.put(None)

    def __init__(self, max_w=sys.maxsize, delay=.01, default_slurm_spec=None):
        self.manager = mp.Manager()
        self.work = self.get_channel()
        self.default_slurm_spec = default_slurm_spec
        self.max_w = max_w
        self.delay = delay
        self.dispatcher = threading.Thread(target=self.dispatch)
        self.dispatcher.start()

    def get_channel(self):
        return self.manager.Queue()

    def close(self):
        self.work.put(None)

    def join(self):
        self.close()
        self.work.get()

    def execute(self, func, args=None, kwargs=None, c=None, slurm_spec=None):
        if args is None: args = tuple()
        if kwargs is None: kwargs = dict()
        if slurm_spec is None: slurm_spec = self.default_slurm_spec
        self.work.put(dill.dumps((slurm_spec, func, args, kwargs, c)))

    @staticmethod
    def exec_wrapper(f, args, kwargs, c, intercom, wid, slurm_spec):
        err = None  # benefit of the doubt
        try:
            if slurm_spec is not None: out = slurm.execute(f, args, kwargs, slurm_spec)
            else: out = f(*args, **kwargs)
        except Exception:
            out, err = None, traceback.format_exc(10)
        if c is not None: c.put((out, err))
        intercom.put(wid)  # notify dispatcher that task is done


if __name__ == '__main__':
    import time

    def log(lc):
        for msg in iter(lc.get,None):
            print(msg)

    def f(x):
        time.sleep((x*3)**.5)
        return str(x**.5)

    def handle_sample(wm, i, lc):
        c = wm.get_channel()
        wm.execute(f, {'x': i}, c=c)
        out, err = c.get()
        if err: raise Exception(err)
        lc.put('%i done (%s)' % (i, str(out)))

    def main():
        slurm_spec = {'cpus-per-task':1,
                      'mem-per-cpu':4000}
        wm = WorkManager(2, 1)
        lc = wm.get_channel()
        lp = mp.Process(target=log, args=(lc,))
        lp.start()
        keep = []
        for i in range(3):
            sh = mp.Process(target=handle_sample, args=(wm, i, lc))
            sh.start()
            keep.append(sh)
        for p in keep: p.join()
        wm.join()
        lc.put(None)
        lp.join()

    main()
