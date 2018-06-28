# Iterate over reads (see a note on input format below), remove ambigous mappings,
# identify and count mismatches to genome and classify reads to converted/unconverted.
#
# If:
# 1) Alignment quality is below min_aq (default=10), or
# 2) CIGAR string is not all M, or
# 3) the number of non T->C errors exceeds max_err_f fraction of bases (default=0.07)
# then read is discarded (into rej_sam).
# Otherwise, read is considered "converted" if there's at least min_cnv_c observed conversions (default=1),
# or "unconverted".
#
# #Ts and #(converted Ts) are collected per read and saved to a 2d histogram file (rt_hist)
# Per every genomic T, the number of reads reporting on it is stored in a wig file (to_wig),
# and the number of times it was observed converted is also written (tc_wig). Also, the
# 2d histogram of observed&genomic genomic Ts is written in gt_hist.
#
# Input is assumed to be a *SORTED* SAM file with original read (=unconverted)
# stored in the query name filed, after a pipe, e.g.:
#  umi:GGGCATTA|ATTGTAGATACTTGACCAATTAGAATTGAATGAGATGAAT   0       chrI    9603    42      40M     *       ...
#
# Variables:
#  Input:
#   genome_tab - a file containing the reference genome with <chr>\t<seq> per line for detecting SNPs
#   debug - "1" to print out read info to stdout, otherwise no printing.
#
#  Parameters: 
#   max_err_f - default % misalignemnt (excluding T->C/A->G) to reject a read. default=0.075.
#   min_aq    - minimal alignment quality for read to not be rejected. default=10;
#   min_cnv_c - the minimal number of converted t->c for the read to be classified as converted
#               default=1.
#  Output:
#   SAM files:
#     rej_sam - a SAM file for rejected reads
#     cnv_sam - a SAM file for converted reads
#     unc_sam - a SAM file for unconverted reads
#     cnv_bed - BED file with unrejected read statistics: chr start end #T:#C 0 strand
#
#   WIG files (all collected only from non-rejected reads):
#     to_wig  - file name (excluding suffix) in which #(observed Ts) per base in the genome
#               is written (2 files are generated <to_wig>.w.wig and <to_wig>.c.wig)
#     tc_wig  - file name (excluding suffix) in which #(converted Ts) per base in the genome
#               is written (2 files are generated <tc_wig>.w.wig and <tc_wig>.c.wig)
#   STATS files:
#     mut_tab - a tab delimited file holdind the complete strand specific mutation matrices
#     rt_hist - a histogram of # of Ts and # converted Ts per read
#     gt_hist - a histogram of # of Ts and # converted Ts per genomic position
#     cl_stats - a tab delimited 3-liner with # reads classified as rejected/converted/unconverted
#

BEGIN {
  # default variable values
  if (genome_tab=="") genome_tab = "/cs/wetlab/genomics/scer/genome/S288C_R64.tab";
  if (debug=="") debug = 0;

  if (max_err_f=="") max_err_f = 0.075;
  if (min_cnv_c=="") min_cnv_c = 1;
  if (min_aq=="") min_aq = 10;

  if (to_wig=="") { to_wig_w = "/dev/null"; to_wig_c = "/dev/null"; }
  else { to_wig_w = to_wig".w.wig"; to_wig_c = to_wig".c.wig"; }
  if (tc_wig=="") { tc_wig_w = "/dev/null"; tc_wig_c = "/dev/null"; }
  else { tc_wig_w = tc_wig".w.wig"; tc_wig_c = tc_wig".c.wig"; }

  if (rt_hist=="") {rt_hist = "/dev/null";}
  if (gt_hist=="") {gt_hist = "/dev/null";}
  if (cl_stats=="") {cl_stats = "/dev/null";}
  if (mut_tab=="") {mut_tab = "/dev/null";}

  if (rej_sam=="") {rej_sam = "/dev/null";}
  if (cnv_sam=="") {cnv_sam = "/dev/null";}
  if (unc_sam=="") {unc_sam = "/dev/null";}
  if (cnv_bed=="") {cnv_bed = "/dev/null";}

  # collect reference genome to memory, future improvement: since input is sorted, can read genomic segments as required
  # future improvement: reverse complement segments once and read directly instead of per read
  while(( getline line < genome_tab) > 0 ) {
    split(line,hs,"\t");
    genome[hs[1]] = hs[2];
  }
  COMP["A"]="T"; COMP["C"]="G"; COMP["G"]="C"; COMP["T"]="A"; COMP["N"]="N";

  rejN = 0;
  cnvN = 0;
  uncN = 0;
}

{ #per row
  split($1,hdr,"|");
  N = split(hdr[2],read,""); #split read
  rpos = $4;
  rchr = $3;
  refseq = substr(genome[rchr],rpos,N);
  is_rev = and(0x10,$2) > 0;
  if (is_rev) {
    tmp = "";
    for (i=N;i>0;i--) {
      tmp = tmp COMP[substr(refseq,i,1)];
    };
    refseq = tmp;
  }

  # count unexpected errors
  err = 0;
  for (i=1;i<=N;i++) {
    ref = substr(refseq,i,1);
    if ( (ref != read[i]) && (ref != "T" || read[i] != "C") ) err++;
  }

  rej = (err >= max_err_f*N) || ($5 < min_aq) || ($6 !~ /^[0-9]+M$/);

  if (debug*1 == 1) {
    print rchr"\t"rpos"\t"is_rev"\t"err"\t"$5"\t"$6"\t"rej"\t"hdr[2]"\t"refseq > "/dev/stderr";
  }

  if ( rej ) {
    print > rej_sam;
    rejN++;
    next; # and stop processing this read further
  }

  r_obsT = 0; r_cnvT = 0;
  for (i=1;i<=N;i++) {
    ref = substr(refseq,i,1);
    if ( ref == "T" ) {
      r_obsT++;
      abspos = is_rev ? rpos+N-i : rpos+i-1; #converting from 1-based indexing to 0 based indexing
      g_obsT[rchr][is_rev][abspos]++; #collect genomic t observed counts
      if ( read[i] == "C" ) {
        r_cnvT++;
        g_cnvT[rchr][is_rev][abspos]++; #collect genomic t conversion counts
      }
    }
    mutmat[is_rev"\t"ref"\t"read[i]]++;
  }
  #print read to bed file
  if ( is_rev ) {
    print rchr"\t"rpos-N"\t"rpos"\t"r_obsT":"r_cnvT"\t"$5"\t-" > cnv_bed;
  }
  else {
    print rchr"\t"rpos"\t"rpos+N"\t"r_obsT":"r_cnvT"\t"$5"\t+" > cnv_bed;
  }

  rhist[r_obsT"\t"r_cnvT]++;
  if ( r_cnvT >= min_cnv_c ) {
    print > cnv_sam; 
    cnvN++;
  }
  else {
    print > unc_sam;
    uncN++;
  }
}

END {
  #print classification statistics 
  print "rejected\t" rejN > cl_stats;
  print "converted\t" cnvN > cl_stats;
  print "unconverted\t" uncN > cl_stats;

  #print mutation statistics
  print "strand\tref\tobs\tcount" > mut_tab;
  for (ref in COMP) for (obs in COMP) print "C\t"ref"\t"obs"\t"mutmat[1"\t"ref"\t"obs]*1 > mut_tab;
  for (ref in COMP) for (obs in COMP) print "W\t"ref"\t"obs"\t"mutmat[0"\t"ref"\t"obs]*1 > mut_tab;

  #print wig files
  for (chr in g_obsT) {
    # watson
    print "variableStep chrom=" chr > to_wig_w;
    print "variableStep chrom=" chr > tc_wig_w;
    if ( 0 in g_obsT[chr] ) {
      for (pos in g_obsT[chr][0]) {
        O = g_obsT[chr][0][pos];
        C = g_cnvT[chr][0][pos]*1;
        print pos" "O > to_wig_w;
        if ( C != 0 ) print pos" "C > tc_wig_w;
        ghist[O"\t"C]++;
      }
    }

    # crick
    print "variableStep chrom=" chr > to_wig_c;
    print "variableStep chrom=" chr > tc_wig_c;
    if ( 1 in g_obsT[chr] ) {
      for (pos in g_obsT[chr][1]) {
        O = g_obsT[chr][1][pos];
        C = g_cnvT[chr][1][pos]*1;
        print pos" "O > to_wig_c;
        if ( C != 0 ) print pos" "C > tc_wig_c;
        ghist[O"\t"C]++;
      }
    }
  }

  #print read 2d histogram
  print "observed\tconverted\tcount" > rt_hist
  for (o_c in rhist) print o_c"\t"rhist[o_c] > rt_hist;

  #print genomic 2d histogram
  print "observed\tconverted\tcount" > gt_hist
  for (o_c in ghist) print o_c"\t"ghist[o_c] > gt_hist;
}
