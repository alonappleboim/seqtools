
#! /bin/bash

# $1 is the list of fastq R1 files to operate on, and is mandatory
# $2 is a file filled with lines of format <bc,sample> where bc is the first 7 bases of R2 (default is "sample_db.csv")
# $3 is the number of reads to operate on in each input file
# $4 is a spike in genome as a bowtie index path, passed to sample handlers
# $5 a stranded TTS bed file - if not "0", a read count is collected per sample, per TTS window (see below)
# $6 TTS upstream distance, default is 300bp
# $7 TTS downstream distance, default is 100bp
# $8 is the script to apply to all samples, default="handle_4tU_sample.bash"

/cs/bd/tools/4tU_scripts/split.bash "$1" $2 $3
/cs/bd/tools/4tU_scripts/run.bash $2 $4 $5 $6 $7 $8
