'''
Written by Alon. June 2017.

Usage examples:

Reads bigwig files from cwd, create hub "DiamideAbcamome2015-test", identifies bigwig files that look like
K3K4me3_TP0_15.bw, sorts by time, groups by Antibody used, and uses coloring instructions from external file:
~/DiamideAbcamome2015_clrs.tab
$ python ~/Dropbox/workspace/seqtools/future/build_hub.py . DiamideAbcamome2015-test "(?P<Ab>[\w\.]+)_TP\d_(?P<time>\d+).bw" -o "time:number" -g "Ab" -f ~/DiamideAbcamome2015_clrs.tab

Reads bigwig files from cwd, create hub "asym-5p_may2017", identifies bigwig files that look like
5324-K36-X-eaf3-m_1.w.bw, groups by strain,res,his,bg,H3, sorts by bg, then H3, then res, then his, then rep, and
finally the strand. Sets default graph type to "bar" in all tracks (no conditions), and adds an indexing entry in the
joint html file.
python ~/Dropbox/workspace/seqtools/future/build_hub.py . asym-5p_may2017 "(?P<strain>[\w]+)(?:-(?P<res>\w+)-(?P<his>\w)-(?P<bg>\w*)-(?P<H3>m|w))?_(?P<rep>\w)_\w+\.(?P<str>w|c)\.bw" -g "strain,res,his,bg,H3" -o "bg:string;H3:string;res:string;his:string;rep:string;str:string" -tp ":graphTypeDefault=bar" -i
'''

import sys
import os
import argparse

# divert to the "right" interpreter TODO: this doesn't work!
# env = '/cs/bd/tools/nflab_env/bin/activate_this.py'
# execfile(env, dict(__file__=env))

URL_BASE = 'http://www.cs.huji.ac.il/labs/nirf/track_hubs'

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input', '-i', type=str, default=None,
                   help='A ucsc track DB file to be converted. Default is stdin.')
    p.add_argument('--output', '-o', type=str, default=None,
                   help='Output file name. ".xml" si added if not present. Default is stdout.')
    p.add_argument('--add_path', '-ap', type=str, default='',
                   help='If given, any url in ucsc hub file that does not start with http (i.e. relative) will be '
                        'added a prefix of %s and this path.')
    args = p.parse_args()
    args.__dict__['input'] = open(args.input) if args.input is not None else sys.stdin
    args.__dict__['output'] = open(args.output,'w') if args.output is not None else sys.stdout
    return args


def parse_track(name, input):
    track = dict(name=name)
    for line in input:
        if not line.strip(): break
        p = line.strip().split(' ')
        track[p[0]] = ' '.join(p[1:])
    if 'aggregate' in track:
        track['subtracks'] = []
    return track


def parse_tracks(input, addpath):
    t, tid = {}, 1
    with input:
        for line in input:
            if line.strip().startswith('track'):
                currt = parse_track(line.strip().split(' ')[1], input)
                if 'color' not in currt: currt['color'] = '0,0,178'
                currt['dispmode'] = 'COLLAPSED'
                if 'graphTypeDefault' in currt and currt['graphTypeDefault'] == 'points':
                    currt['renderer'] = 'LINE_PLOT'
                else:
                    currt['renderer'] = 'BAR_CHART'
                currt['cls'] = 'DataSourceTrack'
                if 'bigDataUrl' not in currt:
                    currt['bigDataUrl'], tid = tid, tid+1
                    currt['cls'] = 'MergedTracks'
                else:
                    if not currt['bigDataUrl'].startswith('http') and addpath:
                        currt['bigDataUrl'] = os.path.sep.join([URL_BASE, addpath, currt['bigDataUrl']])
                if 'parent' in currt:
                    if t['name'] != currt['parent']:
                        names = (currt['name'], t['name'])
                        raise ValueError('track misorder %s inside %s with no parent relationship' % names)
                    t['subtracks'].append(currt)
                else:
                    if t: yield t
                    t = currt


def build_session_xml(tracks):
    yield '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
    yield '<Session genome="sacCer3" hasGeneTrack="true" hasSequenceTrack="true" locus="chrVI:183718-198829" path="" version="8">'
    for line in build_resources(tracks): yield '\t' + line
    for line in build_tracks_panel(tracks): yield '\t' + line
    yield '</Session>'


def build_resources(tracks):
    yield '<Resources>'
    for track in tracks:
        if 'subtracks' in track:
            for st in track['subtracks']:
                yield '\t<Resource path="%s"/>' % st['bigDataUrl']
        else:
            yield '\t<Resource path="%s"/>' % track['bigDataUrl']
    yield '</Resources>'


def build_tracks_panel(tracks):
    yield '<Panel height="496" name="DataPanel" width="1423">'
    for track in tracks:
        for line in build_track(track): yield '\t' + line
    yield '</Panel>'
    yield '<PanelLayout dividerFractions="0.7534039334341907"/>'
    yield '<HiddenAttributes>'
    yield '\t<Attribute name="DATA FILE"/>'
    yield '\t<Attribute name="DATA TYPE"/>'
    yield '\t<Attribute name="NAME"/>'
    yield '</HiddenAttributes>'


def build_track(track):
    str = ('<Track '
           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'altColor="{color}" '
           'autoScale="true" '
           'clazz="org.broad.igv.track.{cls}" '
           'color="{color}" '
           'displayMode="{dispmode}" ' #% COLLAPSED
           'featureVisibilityWindow="-1" '
           'fontSize="10" '
           'id="{bigDataUrl}" '
           'name="{name}" '
           'normalize="false" '
           'renderer="{renderer}" ' #% BAR_CHART
           'sortable="true" '
           'visible="true" '
           'windowFunction="mean" '
           'xsi:type="dataSourceTrack" '
           '>').format(**track)
    if 'subtracks' in track:
        str = str.replace('xsi:type="dataSourceTrack" ', '')
        str = str.replace('autoScale="true" ', 'autoScale="false" ')
        str = str.replace('windowFunction="mean" ', '')
        str = str.replace('normalize="false" ', '')
        str = str.replace('xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ', '')
    yield str
    if 'subtracks' in track:
        for st in track['subtracks']:
            for line in build_track(st): yield '\t' + line
    yield '</Track>'


if __name__ == '__main__':
    args = parse_args()
    tracks = [t for t in parse_tracks(args.input, args.add_path)]
    for line in build_session_xml(tracks): args.output.write(line+'\n')