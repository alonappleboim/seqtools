#! /bin/bash
# The script does two things:
#  1) Generate a TTS window BED file for read counting
#  2) Prints a command per sample to stdout, that can be passed to the parallel_exec.bash script
#
#  It is assumedthat "split.bash" was executed and that the olf (one-line-fasta) files are in the tmp folder.
#
# $1 is a file filled with lines of format <bc,sample> where bc is the first 7 bases of R2 (default is "sample_db.csv")
# $2 is a spike in genome as a bowtie index path, passed to sample handlers
# $3 a stranded TTS bed file - if not "0", a read count is collected per sample, per TTS window (see below)
# $4 TTS upstream distance, default is 200bp
# $5 TTS downstream distance, default is 100bp
# $6 is the script to apply to all samples, default="handle_4tU_sample.bash"

if [ ! -d "tmp" ]; then
  echo "tmp folder does not exist. aborting";
  exit 1;
fi
mkdir STATS 2> /dev/null
mkdir BIGWIG 2> /dev/null
mkdir fastq 2> /dev/null
mkdir BAM 2> /dev/null

db=${1:-sample_db.csv};
spikein=${2:-/cs/wetlab/genomics/klac/bowtie/genome};
TTS=${3:-/cs/wetlab/genomics/scer/annotations/TTS_2018.bed};
TTS_us=${4:-200};
TTS_ds=${5:-100};
script=${6:-/cs/bd/tools/4tU_scripts/handle_4tU_sample.bash};

if [ "$TTS" != "0" ]; then
  ttsw="STATS/tts_window.bed"
  cat $TTS | awk -v ds=$TTS_ds -v us=$TTS_us 'BEGIN {OFS="\t";} {fr=($2-us<0?0:$2-us); to=$2+ds; if ($6=="-") {fr=($2-ds<0?0:$2-ds); to=$2+us}; $2=fr; $3=to; print;}' |\
  sort -k1,1 -k2,2n | cat <(echo "track -$TTS_us $TTS_ds") - > $ttsw
fi

while read entry; do
  IFS=',' read -ra s <<< "$entry"
  bc=${s[0]};
  s=${s[1]};
  echo "passing '$script $s $ttsw $spikein'" > /dev/stderr
  echo "$script $s $ttsw $spikein"
done < "$db"
