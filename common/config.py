INTERPRETER = '/cs/bd/tools/nflab_env/bin/python3.4'

# paths
DATA_PATH = '/cs/bd/SequencingData/'
LOG_PATH = '/cs/bd/logs'
LOCAL_BIN = '/cs/bd/tools/'
URL_BASE = 'http://www.cs.huji.ac.il/labs/nirf/track_hubs'
WWW_PATH = '/cs/bd/track_hubs'

# relative paths
FASTQ_DIR = 'fastq'
BAM_DIR = 'BAM'
BW_DIR = 'BIGWIG'
UNALIGNED_DIR = 'unaligned'
FILTERED_DIR = 'filtered'
TMP_DIR = '.tmp'

# file suffixes
SAM_HDR_SUFF = '.sam.hdr'
BAM_SUFF = '.bam'
TMP_BAM_SUFF = '.tmp.bam'
TMP_BED_SUFF = '.tmp.bed'
TMP_CNT_SUFF = 'cnt.tmp.bed'
FASTQ_SUFF = '.fastq.gz'
UNFILTERED_SUFF = '.unfiltered.bam'

# others
NO_BC_NAME = 'no-barcode'
BC_COUNTS_FNAME = 'bc-counts-%s'

# executables
EXEC = {
    'BOWTIE': 'bowtie2',
    'SAMTOOLS': 'samtools',
    'BEDTOOLS': 'bedtools',
    'SLURM': 'sbatch',
    'BG2W': '/cs/bd/tools/bedGraphToBigWig'
}

#meta
COMMON_GENOMES = {'SCER':
                      {
                          'assembly': 'sacCer3',
                          'bowtie_index': '/cs/wetlab/genomics/scer/bowtie/sacCer3',
                          'chrlens': '/cs/wetlab/genomics/scer/genome/sacCer3.sizes',
                          'annots':
                              {
                                  'tss': '/cs/wetlab/genomics/scer/annotations/weiner2015_tss.tsv',
                                  'tts': '/cs/wetlab/genomics/scer/annotations/weiner2015_tts.tsv'
                              }
                      },
                  'KLAC':
                      {
                          'bowtie_index': '/cs/wetlab/genomics/klac/bowtie/genome',
                          'chrlens': '/cs/wetlab/genomics/klac/genome/xxx',
                          'annots':
                              {
                                  'tss': '/cs/wetlab/genomics/scer/annotations/xx.tsv',
                                  'tts': '/cs/wetlab/genomics/scer/annotations/xx.tsv'
                              }
                       }
                  }

# #error handling
# RETRY_INTERVAL = 5  # sec
# RETRIALS = 3

#meta
SCER_GENOME_LENGTH_PATH = '/cs/wetlab/genomics/scer/genome/sacCer3.sizes'
TTS_MAP = '/cs/wetlab/genomics/scer/annotations/weiner2015_tts.tsv'

STRANDS = {'w': '+', 'c': '-'}

SAMPLEDB_DELIM = ','

def check_third_party():
    """
    check availability of all third party executables listed above with an "_EXEC" suffix
    :return: an error map from an executable
    """
    import subprocess as sp
    import os
    tp = {}
    dn = open(os.devnull)
    for name, val in EXEC.items():
        tp[name] = None
        try: sp.Popen([val, '--help'], stdout=dn, stderr=dn).communicate()
        except OSError as e: tp[name] = e
    return tp


if __name__ == '__main__':
    print(check_third_party())