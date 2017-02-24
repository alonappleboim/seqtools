"""
Given a collection of samples and paths to their BAM files, and an annotation reference file,
generate a corresponding multidimensional matlab file.
"""
from work import *

from common.utils import *
from future.reshape_bam import *


class MainHandler(object):

    def __init__(self, args):
        self.comq = mp.Queue()
        self.a = args
        self.wm = WorkManager(self.comq, self.exec_on=='slurm')

    def split(self):
        """
        generate psuedo-samples by in-silico splitting of data (e.g. by quality, read length, strand,
        etc.)
        :return: None. Updates the samples structure
        """
        pass

    def reduce(self):
        """
        reduce data dimensions - select a specific value in a speicific feature, sum/avg/median/ over
        some feature, etc, ignore certain features, etc.
        :return: None. Updates the samples structure.
        """
        pass

    def hub(self):
        """
        Generate a track hub according to specifications. These include formats and colors, and collating
        of tracks to composite tracks (e.g. overlaying tracks from different time points, or repeats, etc.)
        :return:
        """
        pass

    def reshape_genome(self):
        """
        Reshape genome with respect to certain annotations, according to specifications, specificaly, which
        ReadTransform is used when collecting the reads, and wehter collection is strand-aware or not.
        :return:
        """
    pass


    def export(self):
        """
        Export all data according to specifications.
        :return:
        """
        pass

    def count(self):
        """
        Count reads with respect to certain annotations.
        :return:
        """
    pass

    def meta(self):
        """
        Produce average graphs of sets of annotations according to specifications. These include, locus
        normaliztion (yes/no), which variable is collated on same plot, whether to subplot, and which feature
        is subplotted.
        :return:
        """

    def main(self):
        self.split()
        self.reduce()
        self.hub()
        self.reshape_genome()
        self.count
        token_map = {}
        for s, path in self.split():
            for t in self.transforms:
                self.sample.files['tmp_pkl'] = self.tmp_dir + os.sep + s.short_name() + '.tmp.pkl'
                args = dict(bam_path=s.files['bam'], annot_file=open(self.annot_file),
                            transform=t, out_name='x', output_file=self.sample.files['tmp_pkl'],
                            w=self.w)
                token_map[self.wm.run(reshape, kwargs=args)] = (s, t)
        annots = [line.strip().split('\t')[3] for line in open(self.a.annot_file)]
        dims = (len(annots), self.a.w[1] - self.a.w[0] + 1, len(self.a.transforms)) + \
               tuple([len(f.vals) for f in self.features])
        arr = np.empty(dims)  # a very large array
        leg = OrderedDict()
        for f in self.features: leg[f.name] = sorted(f.vals)
        leg['annot'] = annots
        leg['trans'] = [t.name for t in self.transforms]
        while token_map:
            tok, e, info = self.comq.get()
            sample, transform = token_map.pop(tok)
            if e is not None:
                print('Error while reshaping sample %s with %s '
                      'transform' % (str(sample), transform.name))
                exit()
            else:
                ti = self.transforms.index(transform)
                sidx = tuple(leg[fname].index(sample.fvals[f])
                             for f, fname in zip(self.features, leg.keys()))
                itup = (slice(None), slice(None), ti) + sidx
                arr[itup] = info
        for key, vals in leg.items():
            dtype = np.object if str(vals[0])==vals[0] else np.double
            leg[key] = np.array(vals, dtype=dtype)
        sio.savemat(self.output_file, dict(name=dict(arr=arr, w=self.w, leg=leg)))


def parse_arguments(p):
    args = p.parse_args()
    if args.annot_in is None:
        args.__dict__['annot_in'] = sys.stdin
    else:
        args.__dict__['annot_in'] = open(args.annot_in)
    args.__dict__['w'] = [int(x) for x in args.w[1:-1].split(',')]

    if args.out_name is None:
        args.__dict__['out_name'] = os.path.split(args.output_file)[1].split(os.path.extsep)[0]

    args.__dict__['transform'] = transform_from_string(args.transform)
    return args


def build_parser():
    p = argparse.ArgumentParser()

    g = p.add_argument_group('Input')
    g.add_argument('bam_in', type=str, default=None, help='path to an indexed bam file')
    g.add_argument('--annot_in', '-ain', type=str, default=None,
                   help='path to a "chr,pos,strand" annotation file, default is stdin')
    g.add_argument('--delim', '-d', type=str, default=',', dest='d',
                   help='annotation input delimiter')

    g = p.add_argument_group('Output')
    g.add_argument('--output_file', '-o', type=str, default='reshaped.mat',
                   help='to which output is written.')
    g.add_argument('--out_name', '-on', type=str, default=None,
                   help='name of matrix in matlab struct, default is the suffix of output_name (excluding .mat)')
    g.add_argument('--window', '-w', type=str, default='[-500,250]', dest='w',
                   help='Comma separated limits for the areas around annotations to be collected')
    g.add_argument('--transform', '-T', type=str, default='cov()',
                   help='the transform applied to each read. See -th for details.')

    return p


if __name__ == '__main__':
    c = Collector(parse_arguments(build_parser()))
    c.collect()