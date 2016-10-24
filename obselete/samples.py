import os
import re
import sys
from collections import OrderedDict

from obselete.fastq import merge


class FeatureCollection(OrderedDict):
    def __init__(self):
        super(FeatureCollection, self).__init__()

    def add_feature(self, feat):
        for f in self.values():
            if f.short_name == feat.short_name: raise ValueError()
            if f.name == feat.name: raise ValueError()
        self[feat.name] = feat
        if feat.short_name is None:
            i = 1
            short_names = set(f.short_name for f in self.values())
            while i <= len(feat.name):
                n = feat.name[:i].lower()
                if n in short_names: continue
                feat.short_name = n
                break
        if feat.short_name is None:
            raise ValueError()


class Feature(object):
    def __init__(self, name, type, short_name=None, units=None):
        self.name = name.lower()
        self.short_name = short_name
        self.strtype = type.lower()
        self.type = str if type == 'str' else int if type == 'int' else float
        self.units = units
        self.vals = set([])

    def __str__(self):
        return '%s(%s)[%s]:%s' % (self.name, self.short_name, self.units, self.strtype)

    def __repr__(self): return str(self)


class Sample(object):
    def __init__(self, context):
        self.fvals = OrderedDict()
        self.context = context
        self.barcode = None
        self.files = {}

    def fastq(self):
        return context.fastq_dir + os.sep + self.base_name() + '.fastq.gz'

    def base_name(self):
        return '_'.join('%s-%s' % (f.short_name, str(v)) for f, v in self.fvals.items())

    def full_name(self):
        return '_'.join('%s-%s' % (f.name, str(v)) for f, v in self.fvals.items())

    def

    def handle(self, in1, in2, context):
        if context.start_from <= FASTQ:
            fastq = self.fastq()
            try:
                kwargs = dict(pre1=in1, pre2=in2, fastq=fastq, umi_len=context.umi_len, bc_len=context.bc_len)
                n = context.w_manager.execute(merge, kwargs=kwargs)
                context.statq.put((self, n))
                context.statusq.put((self, FASTQ, None))
            except e:
                context.statusq.put((self, FASTQ, e))





    def __repr__(self):
        return self.base_name()

    def __hash__(self):
        return hash(tuple(self.fvals.values()))


class SampleDBParser(object):
    FEAT_PAT = re.compile('\s*(?P<name>\w+)\s*(?:\((?P<short_name>\w+)\))?'
                          '\s*:(?P<type>\w+)(:?\[(?P<units>\w+)\])?')

    def __init__(self, delim=',', n_samples=sys.maxsize):
        self.d = delim
        self.n = n_samples
        self.err = []
        self.info = []

    def _parse_header(self, hdr):
        hdr = hdr.split(self.d)
        if hdr[0] != 'barcode':
            self.err.append('first column in sample db needs to be the "barcode" column')
        features = FeatureCollection()
        f_pos_map = {}
        for i, f in enumerate(hdr[1:]):
            if f.startswith('#'):
                self.info.append("ignoring column %s in sample db" % f)
            else:
                m = re.match(SampleDBParser.FEAT_PAT, f)
                if m is None:
                    self.err.append("couldn't understand feature '%s' in sample_db file, format should be: "
                                    "<name>(<short_name>):(str|int|float)[units] (short_name and units are "
                                    "optional) or column is ignored if it starts with '#'" % f)
                else:
                    try:
                        f_pos_map[i] = Feature(**m.groupdict())
                    except ValueError:
                        snames = '\n'.join(f.short_name for f in features.values)
                        self.err.append("features must have distinct names and short_names - %s "
                                        "appears at least twice (or its short_name matched a "
                                        "previous generated short name):\n%s" % f, snames)
                    features.add_feature(f_pos_map[i])
        return features, f_pos_map

    def _parse_samples(self, db, fmap):
        b2s, bc_len = OrderedDict(), None
        for i, line in enumerate(db):
            if line.strip()[0] == '#': continue  # comment
            if i >= self.n: break  # limit number of samples
            sample = Sample()
            for j, val in enumerate(line.strip().split(self.d)):
                val = val.strip()
                if j == 0:
                    if bc_len is None: bc_len = len(val)  # first barcode
                    if bc_len != len(val):
                        self.err.append("barcode %s has a different length" % val)
                    if val in b2s:
                        self.err.append("barcode %s is not unique" % val)
                    sample.barcode = val
                elif j in fmap:
                    f = fmap[j]
                    try:
                        v = f.type(val)
                        f.vals.add(v)
                        sample.fvals[f] = v
                    except ValueError:
                        self.err.append("couldn't cast value %s in sample %i, feature '%s' to "
                                        "given type - %s." % (val, i + 1, f.name, f.strtype))
            if hash(sample) in [hash(s) for s in b2s.values()]:
                self.err.append("2 samples (or more) seem to be identical - %s" % sample)
            b2s[sample.barcode] = sample
        return b2s, bc_len

    def parse(self, path):
        self.errs.clear()
        self.info.clear()
        with open(path) as sdb:
            exp = re.match('.*experiment.*:\s+(\w+)', sdb.readline())
            if exp is None:
                raise ValueError('barcodes file should contain a header with experiment name: '
                                 'experiment: <expname>')
            features, f_pos_map = self._parse_header(sdb.readline())
            if self.err: raise ValueError('\n'.join(self.err))
            b2s, bc_len = self.parse_samples(sdb, f_pos_map)
            if self.err: raise ValueError('\n'.join(self.err))
        return bc_len, b2s, features


