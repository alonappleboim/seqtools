#! /bin/awk
#  prepare reads in fastq format for alignment to a converted (ACG-only) genome
#
#  Variables:
#    strip_a - strip a's or not strip a
#    min_rlen - minimal remaining read length after stripping for read to be printed. default=25.
#    stripc_file - number of rejected reads is written to this file, if given
#    ahist_file - histgram of #trailing As is written to this file, if given
#
BEGIN {
  # default variable values
  if (strip_a == "") strip_a = 0;
  if (min_rlen == "") min_rlen = 25;
  if (stripc_file == "") stripc_file = "/dev/null";
  if (ahist_file == "") ahist_file = "/dev/null";
  rej = 0;
}
{
  s = $2;
  q = $4;
  tmp = s;
  gsub(/A*$/,"",tmp);
  ah[length(s)-length(tmp)]++; # count #trailing As into histogram ah
  if (strip_a == 1) s = tmp;
  if (length(s) < min_rlen) { rej++; next; } # count as rejected and do not print the read
  q = substr($4,1,length(s));
  w = s; gsub(/T/,"C",w);
  print $1"|"s"\n"w"\n+\n"q;
}

END {
  print rej > stripc_file;
  for (i in ah) print i"\t"ah[i] > ahist_file;
}
