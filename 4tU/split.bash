#! /bin/bash

# $1 is the list of fastq R1 files to operate on, and is mandatory
# $2 is a file filled with lines of format <bc,sample> where bc is the first 7 bases of R2 (default is "sample_db.csv")
# $3 is the number of reads to operate on in each input file

mkdir tmp 2> /dev/null

SPL="/cs/bd/tools/scripts/split_barcodes.py"

db=${2:-sample_db.csv};
nr=${3:-1000000000};
echo "splitting barcodes..." > /dev/stderr
for f1 in $1; do
  f2=${f1/R1/R2};
  pref=`basename "$f1"`
  pref=${pref%%.*}
  pref=${pref##_}
  echo $pref $f1 $f2 > /dev/stderr
  echo "  in file $f1, outputing to tmp/hd*$pref" > /dev/stderr
  paste <(zcat $f1) <(zcat $f2) | paste - - - - | head "-$nr" | \
  python $SPL $db -d 0 -od "tmp/hd0$pref" |\
  python $SPL $db -d 1 -od "tmp/hd1$pref" > /dev/null
done
