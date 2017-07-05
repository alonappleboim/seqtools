import sys
import argparse
import os
import subprocess as sp
import shutil
import re
from collections import OrderedDict
from matplotlib.pyplot import cm
import numpy as np

DEFAULT_CMAP = 'rainbow'

# divert to the "right" interpreter TODO: this doesn't work!
env = '/cs/bd/tools/nflab_env/bin/activate_this.py'
execfile(env, dict(__file__=env))

TYPES = {'string': str, 'number': float}

URL_BASE = 'http://www.cs.huji.ac.il/labs/nirf/track_hubs'
TRACKS_CENTRAL = '/cs/bd/track_hubs'
INDEX_FILE = '/cs/bd/track_hubs/index.html'


def shades(clr, N):
    for i in range(N):
        yield [int((255 * i / N + gc) / 2) for gc in clr]


class TrackProps(OrderedDict):

    def __init__(self, *args, **kwargs):
        super(TrackProps, self).__init__(*args, **kwargs)

    def format(self, props=[], pref=''):
        p = self.copy()
        for k,v in props: p[k] = v
        return '\n'.join('%s%s %s' % (pref, k, v) for k,v in p.items())


DEF_PARENT_PROPS = TrackProps([('track','none'),
                               ('container','multiWig'),
                               ('aggregate','transparentOverlay'),
                               ('type','bigWig'),
                               ('autoScale','on'),
                               ('visibility','full'),
                               ('shortLabel','none'),
                               ('longLabel','none'),
                               ('maxHeightPixels','100:32:8'),
                               ('priority','50')])


DEF_TRACK_PROPS = TrackProps([('track','none'),
                              ('parent','none'),
                              ('type','bigWig 0 1000'),
                              ('color','255,255,255'),
                              ('altColor','255,255,255'),
                              ('alwaysZero','on'),
                              ('yLineOnOff','on'),
                              ('visibility','full'),
                              ('smoothingWindow','4'),
                              ('windowingFunction','mean'),
                              ('graphTypeDefault','points'),
                              ('bigDataUrl','none')])


def get_cmap(cname, n):
    cmap = []
    if os.path.exists(cname):
        try:
            for i, line in enumerate(open(cname)):
                cmap.append([int(x) for x in line.strip().split('\t')])
                if i == n: break
            if i < n:
                raise ValueError('Not enough colors, reverting to default colors')
        except Exception as e:
            sys.stderr.write('Could not parse RGB file: %s, %s, reverting to default colors' % (cname, str(e)))
            cmap = []
            cname = DEFAULT_CMAP
    if not cmap:
        try:
            cmapobj = cm.get_cmap(cname)
            cmap = np.asarray(cmapobj(np.linspace(0, 1, n))[:, :3]*255, dtype=int)
        except Exception as e:
            raise ValueError(sys.stderr.write('Could not get colormap %s (%s)' % (cname, str(e))))
    return cmap


def pprint_rgb(c):
    return '('+','.join(str(x) for x in c) + ') #' + ''.join(hex(int(x))[2:] for x in c)


class Track(object):

    def __init__(self, src, vars):
        self.src = src
        self.vars = vars
        self.color = None
        self.aprops = OrderedDict()

    def __repr__(self):
        return ','.join(['%s-%s' %(k,v) for k,v in self.vars.items()])+ ' #' + ''.join(hex(x)[2:] for x in self.color)

    def set_props(self, condprops):
        for conds, props in condprops:
            passed = True
            for var, val in conds:
                if self.vars[var] != val:
                    passed = False
                    break
            if passed:
                for k,v in props: self.aprops[k] = v


    def get_name(self):
        return '_'.join(['%s-%s' % (k, v) for k, v in self.vars.items()])


class Hub(object):

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.gtracks = OrderedDict()
        self.vars = OrderedDict()
        self.grouped_by = []
        self.var_order = []

    def organize_by_goc(self, goc_path):
        clrs = []
        with open(goc_path) as GOC:
            vars = GOC.readline().strip.split('\t')
            gtracks = OrderedDict()
            assert (vars[-1] == 'RGB', 'Last column in GOC file should be "RGB"')
            vars = vars[:-1]
            for v in vars: assert(v in self.vars)
            for line in GOC:
                tg = []  # track group
                vals = line.strip().split('\t')
                vals, c = vals[:-1], vals[-1]
                for ts in self.gtracks.values():
                    for t in ts:
                        add = True
                        for var, val in zip(vars, vals):
                            if t.vars[var] != val:
                                add = False
                                break
                        if add: tg.append(t)


    def from_path(self, path, ve):
        path = os.path.abspath(path)
        for f in os.listdir(path):
            m = ve.match(f)
            if not m: continue
            m = tuple(m.groupdict().items())
            for var, val in m:
                if var not in self.vars: self.vars[var] = set([])
                self.vars[var].add(val)
            self.gtracks[m] = [Track(src=path+os.path.sep+f, vars=dict(m))] # single track groups to begin with
        self.grouped_by = list(self.vars.keys())

    def update_track_props(self, cond_props):
        for ts in self.gtracks.values():
            for t in ts: t.set_props(cond_props)

    def regroup_tracks(self, group_by):
        if not group_by: return
        gtracks = {}
        for ts in self.gtracks.values():
            for t in ts:
                tvars = dict(t.vars)
                k = tuple([(gb, tvars[gb]) for gb in group_by])
                if k not in gtracks: gtracks[k] = []
                gtracks[k].append(t)
        self.gtracks = gtracks
        self.grouped_by = group_by

    def sort_tracks(self, order_by):
        ordered = self.gtracks.copy()
        for var, typ in order_by[::-1]:
            self.var_order = [var] + self.var_order
            self.vars[var] = sorted(self.vars[var], key=lambda x:TYPES[typ](x))
            if var not in self.vars: raise ValueError('Unknown variable: %s' % var)
            if var not in self.grouped_by: # inner group sorting
                ordered = OrderedDict((gvars, sorted(tracks, key=lambda x: TYPES[typ](x.vars[var])))
                                      for gvars, tracks in ordered.items())
            else:  # group sorting
                ordered = OrderedDict((gvars, ordered[gvars])
                                      for gvars in sorted(ordered.keys(), key=lambda x:TYPES[typ](dict(x)[var])))
        self.gtracks = ordered

    def color_tracks(self, cmap_name, color_by):
        if not color_by: color_by = list(self.grouped_by)
        for c in color_by: assert(c in self.grouped_by) # only outer-grouping variables can be in color_by
        ckeys, gkeys = {}, []
        for gvars, tracks in self.gtracks.items():
            gvars = dict(gvars)
            key = tuple([gvars[var] for var in color_by])
            gkeys.append(key)
            if key in ckeys: continue
            ckeys[key] = len(ckeys)
        cmap = get_cmap(cmap_name, len(ckeys))
        for key, (gvars, tracks) in zip(gkeys, self.gtracks.items()):
            gclr = cmap[ckeys[key]][:3]
            for track, clr in zip(tracks, shades(gclr,len(tracks))): track.color = clr

    def deploy(self, dst):
        dstpath = os.path.abspath(TRACKS_CENTRAL + os.path.sep + dst)
        sys.stderr.write('Generating hub at %s...\n' % dstpath)
        if not os.path.isdir(dstpath): os.mkdir(dstpath)
        else:
            raise IOError('destination exists, aborting.')
            # delete = input('destination folder exists, delete and overwrite? (N/y) ')
            # if delete != 'y': raise IOError('destination exists, aborting.')
            # shutil.rmtree(dstpath)
            # os.mkdir(dstpath)
        with open(dstpath + os.path.sep + 'hub.txt', 'w') as HUB:
            HUB.write('\n'.join(["hub %s" % self.name,
                                 "shortLabel %s" % self.name,
                                 "longLabel %s" % self.full_name,
                                 "genomesFile genomes.txt",
                                 "email %s" % self.email]))
        with open(dstpath + os.path.sep + 'genomes.txt', 'w') as GEN:
            GEN.write('\n'.join(["genome %s" % args.genome_assembly,
                                 "trackDb trackDB.txt"]))
        with open(dstpath + os.path.sep + 'trackDB.txt', 'w') as T:
            for gvars, tracks in self.gtracks.items():
                gvars = dict(gvars)
                gname = '_'.join(['%s-%s' % (v, gvars[v]) for v in self.grouped_by])
                sys.stderr.write('Handling %s\n' % gname)
                hdr = DEF_PARENT_PROPS.format([('shortLabel',gname),('longLabel',gname),('track',gname)])
                T.write(hdr + '\n\n')
                for t in tracks:
                    tname = t.get_name()
                    tdst = dstpath + os.path.sep + tname + '.bw'
                    shutil.copy(t.src, tdst)
                    if args.link:
                        pass #remove old file and replace with a link to the new version in the track hub
                    turl = os.path.sep.join([URL_BASE, dst, tname + '.bw'])
                    cstr = ','.join(str(i) for i in t.color)
                    props = [('parent',gname), ('track',tname), ('color',cstr), ('altColor',cstr), ('bigDataUrl',turl)]
                    props = list(t.aprops.items()) + props
                    tentry = DEF_TRACK_PROPS.format(props, '\t')
                    T.write(tentry + '\n\n')
        sp.call('chmod -R 777 %s' % dstpath, shell=True)
        return os.path.sep.join([URL_BASE, dst, 'hub.txt'])

def parse_cond_props(condpropstr):
    # parse, e.g.: strain=BY+time=0:type=bars+color=10,40,190;strain=mut:visibility=hide'
    cps = []
    if condpropstr:
        for condprop in condpropstr.split(';'):
            conds, props = condprop.split(':')
            conds = [c.split('=') for c in conds.split('+')] if conds else []
            props = [c.split('=') for c in props.split('+')] if props else []
            cps.append((conds, props))
    return cps


def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument('bw_path', type=str, help='BigWig folder to hubify')
    p.add_argument('name', type=str, help='Hub name')
    p.add_argument('--output', type=str, help='destination folder name, default is "name"', default=None)
    p.add_argument('--full_name', '-fn', type=str, help='Hub full name', default=None)
    p.add_argument('--email', '-e', type=str, help='contact email', default='address@email.domain')
    p.add_argument('var_regexp', type=str,
                   help=('Name matching python regexp that also defines the variables. '
                         'e.g. (?P<mod>\w+)_(?P<time>\d+)\.bw'))
    p.add_argument('--group_by', '-g', type=str, default=None,
                   help=('Variarble name(s) by which track are grouped, i.e. shown on the same track, '
                         'different shades of the same color'))
    p.add_argument('--order_by', '-o', type=str, default=None,
                   help=('A list of variarbles and how to sort them. order of variables determines order in resulting'
                         ' hub. e.g. "mod:str;time:int" will be first sorted by "mod"(string), then by "time"(number)'))
    p.add_argument('--color_by', '-c', type=str, default=None,
                   help=('A comma-separated list by which colors are allocated. default is all variables'))
    p.add_argument('--colors', type=str, default=DEFAULT_CMAP,
                   help=('A string corresponding to a matplotlib colormap'))
    p.add_argument('--goc_file', '-f', type=str, default=None,
                   help=('A path to a tab-delimited file with variable columns, and an RGB column. The variables are '
                         'the grouping variables, the order is used to order the groups. Tracks with variable-values '
                         'that are not present in the file are filtered out, and colors will be the group colors.'
                         'Replaces and overrides the group_by/order_by/color_by/colors options (not that inner-group '
                         'order is still used'))
    p.add_argument('--link', '-l', action='store_true',
                   help=('Whether original tracks should be replaced by sof links to newly copied files')) #TODO
    p.add_argument('--track_props', '-tp', type=str, default='',
                   help=('By default, all tracks are generated with these properties: \n%s\n'
                         'To override any property or add properties to track groups, give a '
                         '"conditions:props;conditions:props.." clause. by default props are applied to all tracks, '
                         'unless they did not comply with condition. Conditions are handled in order, so if some '
                         'property of a track is changed by several conditions the last condition wins. For example:\n'
                         '"strain=BY+time=0:type=bars+color=10,40,190;strain=mut:visibility=hide') % DEF_TRACK_PROPS.format())
    p.add_argument('--not_verbose', '-nv', action='store_true',
                   help=('Whether the hub description file should not be printed to standard output.'))
    p.add_argument('--genome_assembly', '-ga', type=str, default='sacCer3',
                   help=('The standard genomic assembly name in UCSC browser'))
    p.add_argument('--index', '-i', action='store_true',
                   help=('Add hub to index.html file, you will prompted for further details'))
    args = p.parse_args()
    if args.full_name is None: args.__dict__['full_name'] = args.name
    args.__dict__['group_by'] = args.group_by.split(',') if args.group_by is not None else []
    args.__dict__['order_by'] = [x.split(':') for x in args.order_by.split(';')] if args.order_by is not None else []
    args.__dict__['color_by'] = args.color_by.split(',') if args.color_by is not None else []
    args.__dict__['var_regexp'] = re.compile(args.var_regexp)
    args.__dict__['track_props'] = parse_cond_props(args.track_props)
    if args.output is None: args.__dict__['output'] = args.name
    return args


def add_to_index(hubname, url):
    desc = raw_input('Please enter a hub description (html compatible, e.g. urls can be added as <a href="url">"link name"</a>):\n')
    entry = ['  <tr>',
             '    <td> <a href="%s">%s</a></td>' % (url, hubname),
             '    <td>%s</td>' % desc,
             '  <td>%s</td>' % hubname,
             '  </tr>']
    original = []
    with open(INDEX_FILE) as F:
        last_tr = -1
        for i, line in enumerate(F):
            original.append(line.strip())
            if '<tr>' in line: last_tr = i
    updated = original[:last_tr-1] + entry + original[last_tr:]
    with open(INDEX_FILE,'w') as F:
        F.write('\n'.join(updated))
    sys.stderr.write('updated %s!' % INDEX_FILE)


if __name__ == '__main__':
    args = parse_args()
    hub = Hub(name=args.name, full_name=args.full_name, email=args.email)
    hub.from_path(args.bw_path, args.var_regexp)
    if args.goc_file is not None:
        hub.organize_by_goc(args.goc_file)
        hub.sort_tracks(args.order_by)
    else:
        hub.regroup_tracks(args.group_by)
        hub.sort_tracks(args.order_by)
        hub.color_tracks(args.colors, args.color_by)
    hub.update_track_props(args.track_props)
    url = hub.deploy(args.output)
    sys.stderr.write('Hub available at %s\n' % url)
    if args.index:
        add_to_index(args.name, url)