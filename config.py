
# paths
DATA_PATH = '/cs/bd/SequencingData/TranSEQ/'
LOG_PATH = '/cs/bd/TranSEQ/logs'
LOCAL_BIN = '/cs/bd/tools/'

# output
UNALIGNED_NAME = 'unaligned'
TMP_NAME = '.tmp'
NO_BC_NAME = 'no-barcode'
BC_COUNTS_FNAME = 'bc-counts-%s'
SAM_HDR_SUFF = '.sam.hdr'
TMP_BAM_SUFF = '.tmp.bam'
BT_STATS_SUFF = '.bowtie.stats'

# executables
BOWTIE_EXEC = 'bowtie2'
SAMTOOLS_EXEC = 'samtools'
BG2W_EXEC = LOCAL_BIN+'bedGraphToBigWig'
GC_EXEC = 'genomeCoverageBed'
