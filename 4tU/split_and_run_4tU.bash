#! /bin/bash

# $1 is the list of fastq R1 files to operate on, and is mandatory
# $2 is a file filled with lines of format <bc,sample> where bc is the first 7 bases of R2 (default is "sample_db.csv")
# $3 is the number of reads to operate on in each input file
# $4 is a spike in genome as a bowtie index path, passed to sample handlers
# $5 a stranded TTS bed file - if not "0", a read count is collected per sample, per TTS window (see below)
# $6 TTS upstream distance, default is 300bp
# $7 TTS downstream distance, default is 100bp
# $8 is whether splitting part should be skipped (i.e. tmp exists with split reads)
# $9 is the script to apply to all samples, default="handle_sample.bash"

mkdir STATS 2> /dev/null
mkdir BIGWIG 2> /dev/null
mkdir fastq 2> /dev/null
mkdir tmp 2> /dev/null
mkdir BAM 2> /dev/null

SPL="/cs/bd/tools/scripts/split_barcodes.py"

db=${2:-sample_db.csv};
nr=${3:-1000000000};
spikein=${4:-/cs/wetlab/genomics/klac/bowtie/genome};
TTS=${5:-/cs/wetlab/genomics/scer/annotations/TTS_2018.bed};
TTS_us=${6:-300};
TTS_ds=${7:-100};
skip_split=${8:-0};
comm=${9:-/cs/bd/tools/4tU_scripts/handle_4tU_sample.bash};

if [ "$skip_split" -ne "1" ]; then
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
fi

echo $TTS > /dev/stderr

if [ "$TTS" != "0" ]; then
  ttsw="STATS/tts_w.bed"
  cat $TTS | awk -v ds=$TTS_ds -v us=$TTS_us 'BEGIN {OFS="\t"; print "track -"us" "ds;} {fr=($2-us<0?0:$2-us); to=$2+ds; if ($5=="-") {fr=($2-ds<0?0:$2-ds); to=$2+us}; $2=fr; $3=to; print;}' > $ttsw
fi

echo "handling samples..." > /dev/stderr
while read entry; do
  IFS=',' read -ra s <<< "$entry"
  bc=${s[0]};
  s=${s[1]};
  echo "  passing '$comm $s $ttsw $spikein'" > /dev/stderr
  echo "$comm $s $ttsw $spikein"
done < "$db"
