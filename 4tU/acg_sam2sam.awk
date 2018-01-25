#! /bin/awk
#  convert a three-letter-genome SAM to a regular SAM
#  Note that unaligned reads are discarded!
#
#  Variables:
#    chr_len - a tab delimited mapping of chr names to lengths
#
BEGIN {
  # default variable values
  if (genome_lengths=="") genome_lengths = "/cs/wetlab/genomics/scer/genome/sacCer3.sizes"
  # collect reference genome lengths to memory
  while(( getline line < genome_lengths) > 0 ) {
    split(line,hs,"\t");
    gl[hs[1]] = hs[2];
  }
  OFS = "\t";
}
{
# umi:CGAGGTTG|ATTCAGTAATGAGCACAAAAGTTGACATTCATCCCATTTT   0       chrI-w  14010   37      40M     *       0       0       ACCCAGCAACGAGCACAAAAGCCGACACCCACCCCACCCC
  if ($3 == "*") { print; next; } #unaligned, nothing to do
  split($3,tmp,"-")
  c = tmp[1];
  s = tmp[2];
  $3 = c; #update chromosome to original chromsome name
  if (s=="c") { 
    #if crick, update flag and position
    $2 += 16
    p = gl[c]-($4-1)-(length($10)-1); #len(c)-(len(read)-1)-(start-1)
    $4 = p;
#    print $4;
  }
  print; 
}
