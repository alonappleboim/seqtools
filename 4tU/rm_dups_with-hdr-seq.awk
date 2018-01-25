# remove duplicates based on their umi, start position, and alignment cigar string, assuming input is sorted
#  dstat_file = file to which detailed duplication statistics are written
#  dhist_file = file to which duplicate count histogram is written  

# HANDLE THE CASE THAT in which the read header contains the original sequence  (umi|<seq>)

BEGIN {
    if (dstat_file=="") { dstat_file = "/dev/null"; }
    if (dhist_file=="") { dhist_file = "/dev/null"; }
}

{
    umi = substr($1,1,13);
    if (prev != $4) {
        for (x in u) print u[x];
        for (x in c) {
            print x,c[x] > dstat_file;
            h[c[x]]++;
        }
	delete u;
	delete c;
	c[umi"\t"$4"\t"$6] ++;
	u[umi"\t"$4"\t"$6] =$0; #$4 = alignment position, $9 fragment length
    }
    else {
        c[umi"\t"$4"\t"$6] ++;
	u[umi"\t"$4"\t"$6] = $0;
    }
    prev = $4;
}
END {
        for (x in u) print u[x];
        for (x in c) print x,c[x] > dstat_file; 
        for (y in h) print y"\t"h[y] > dhist_file;
} #print last unique read




