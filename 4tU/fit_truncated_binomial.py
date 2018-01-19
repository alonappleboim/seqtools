import argparse
import sys
import pickle

from scipy.stats import binom
from scipy import optimize
import numpy as np


def fit_truncated_binomial(input, ts=np.arange(5), poiss_noise=0.01, N=100):
    hist = np.loadtxt(input, skiprows=1, dtype=float)
    phat = np.empty((len(ts), N))
    phat[:] = np.nan
    for ti, t in enumerate(ts):
        dtf = hist[hist[:, 1] >= t, :]  # truncate
        n, x = dtf[:, 0], dtf[:, 1]
        for ni in range(N):
            c = np.random.poisson(dtf[:, 2].T + poiss_noise)
            # negative log likelihood of the t-truncated binomial distribution:
            ttbino_nll = lambda p: -np.dot(c, binom.logpmf(x, n, p) - np.log(1 - binom.cdf(t - 1, n, p)))
            phat[ti, ni] = optimize.fminbound(ttbino_nll, 0, 1)
    return phat, dict(ts=ts, N=N, poiss_noise=poiss_noise)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('input', type=str,
                   help='path to a sparse 2D histogram counting file with one header line and 3 tab delimited values:'
                        '1) #observed Ts 2) #converted Ts 3) #reads')
    p.add_argument('--truncate_ts', '-tt', type=str, help='a list of comma-separated integers to truncate and fit',
                   default='0,1,2,3,4', dest='ts')
    p.add_argument('--bootstrap_n', '-N', type=int, help='number of poisson-noised bootstrap repetitions to perform',
                   default=100, dest='N')
    p.add_argument('--poiss_noise', '-pn', type=float, help='poisson-noise added to count data for bootstrapping',
                   default=.01)
    p.add_argument('--output', '-o', type=str, default=None, help='output file name, default is stdout')
    args = p.parse_args()

    args.__dict__['output'] =  open(args.output, 'wb') if args.output else sys.stdout
    args.__dict__['ts'] = [int(x) for x in args.ts.split(',')]

    return args


if __name__ == '__main__':
    args = parse_args()
    phat, params = fit_truncated_binomial(args.input, ts=args.ts, poiss_noise=args.poiss_noise, N=args.N)
    hdr = pickle.dumps(params)
    np.savetxt(args.output,phat,header=hdr)

