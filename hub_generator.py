
def write_hub(pAs, hdr, hub_path):
    print 'writing to: ' + hub_path
    if not os.path.exists(hub_path): os.mkdir(hub_path)
    paths = [hub_path+os.path.sep+x for x in hdr]
    handles = [(open(p+'_F.bg', 'wb'),open(p+'_R.bg', 'wb')) for p in paths]
    # 	for h, hd in zip(handl es,hdr):
    # 		h.write('track type=bedGraph name="Wilkening2013 - %s" description="poly A sites - %s" '
    # 				'visibility=full color=200,100,0 altColor=0,100,200 priority=20\n' %  (hd, hd))
    # collect data
    for i, (k,v) in enumerate(pAs.items()):
        if i % 5e4 == 0: print 'At %i...' % i
        chr, start, end, strand = k
        chr = 'chr' + (toRoman(int(chr)) if chr != '17' else 'M')
        end = str(int(end)+1)
        for j, h in enumerate(handles):
            if v[j] == '0': continue
            val, hi = ('-'+v[j], 1) if strand == '0' else (v[j], 0)
            h[hi].write(' '.join([chr, start, end, val])+'\n')

    sacpath = hub_path + os.path.sep + 'sacCer3'
    if not os.path.exists(sacpath): os.mkdir(sacpath)
    trackfile = open(sacpath + os.path.sep + 'trackDB.txt', 'wb')
    #convert to bigwig and write into trackDB
    for (hf, hr), p, hd in zip(handles, paths, hdr):
        print 'Sorting %s...' % hd
        hf.close()
        hr.close()
        sp.call('sort -k1,1 -k2,2n %s > %s_sorted.bg' % (p+'_F.bg',p+'_F'), shell=True)
        sp.call('sort -k1,1 -k2,2n %s > %s_sorted.bg' % (p+'_R.bg',p+'_R'), shell=True)
        print 'And writing...'
        sp.call('~/Dropbox/software/%s/bedGraphToBigWig %s_F_sorted.bg '
                '%s %s_F.bw' % (bindir, p,CHRLEN_PATH,p), shell=True)
        sp.call('~/Dropbox/software/%s/bedGraphToBigWig %s_R_sorted.bg '
                '%s %s_R.bw' % (bindir, p, CHRLEN_PATH, p), shell=True)
        os.remove(p+'_F.bg')
        os.remove(p+'_R.bg')
        os.remove(p+'_F_sorted.bg')
        os.remove(p+'_R_sorted.bg')
        entries = (hd+'_F', URL+hd+'_F.bw', 'pA - %s (F)'%hd, 'Wilkening 2013 polyA data, %s (forward)' % hd)
        trackfile.write("track %s\n"
                        "bigDataUrl %s\n"
                        "shortLabel %s\n"
                        "longLabel %s\n"
                        "type bigWig\n"
                        "visibility full\n"
                        "viewLimits 0:500\n\n" % entries)
        entries = (hd+'_R', URL+hd+'_R.bw', 'pA - %s (R)'%hd, 'Wilkening 2013 polyA data, %s (reverse)' % hd)
        trackfile.write("track %s\n"
                        "bigDataUrl %s\n"
                        "shortLabel %s\n"
                        "longLabel %s\n"
                        "type bigWig\n"
                        "visibility full\n"
                        "viewLimits -500:0\n\n" % entries)
    print 'Setting up hub...'
    hubfile = open(hub_path+os.path.sep+'hub.txt','wb')
    hubfile.write("hub %s\n"
                  "shortLabel pA sites\n"
                  "longLabel Data relevant to 3' processing and alternative UTRs\n"
                  "genomesFile genomes.txt\n"
                  "email  alonappleboim@gmail.com" % os.path.split(hub_path)[-1])
    genomesfile = open(hub_path+os.path.sep+'genomes.txt','wb')
    genomesfile.write("genome sacCer3\n"
                      "trackDb sacCer3%strackDB.txt" % os.path.sep)