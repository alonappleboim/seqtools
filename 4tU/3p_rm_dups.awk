# remove duplicates based on their umi, start position, and alignment cigar string, assuming input is sorted
# vars:
#  dstat_file = file to which detailed duplication statistics are written
#  dhist_file = file to which duplicate count histogram is written  
#
# for example:
# bowtie2 -1 $f1 -2 $f2 -x $index_path| samtools sort | awk -v dcnt_file=dup.cnts -f /cs/bd/tools/scripts/rm_dups.awk| samtools view -t $path_to_chr_lengths -b > $bam_file
# samtools index $bam_file
# cat dup.cnts | awk '{print $3}' | sort | uniq -c > dup.hist
#

BEGIN {
    if (dstat_file=="") {dstat_file = "/dev/null";}
    if (dhist_file=="") {dhist_file = "/dev/null";}}
{
    if (prev != $4) {
        for (x in u) print u[x];
	for (x in c) {
            print x,c[x] > dstat_file;
            h[c[x]]++;
        }
	delete u;
	delete c;
	c[$1"\t"$4"\t"$6]++;
	u[$1"\t"$4"\t"$6] =$0; #$4 = alignment position, $9 fragment length
    }
    else {
        c[$1"\t"$4"\t"$6]++;
	u[$1"\t"$4"\t"$6] = $0;
    }
    prev = $4;
}
END {
	for (x in u) print u[x];
	for (x in c) print x,c[x] > dstat_file; 
        for (y in h) print y"\t"h[y] > dhist_file;
} #print last unique read

