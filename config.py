import sys

INTERPRETER = '/usr/local/opt/python3/bin/python3.5'
if not sys.executable == INTERPRETER:  # divert to the "right" interpreter
    import subprocess as sp
    import os
    scriptpath = os.path.abspath(sys.modules[__name__].__file__)
    sp.Popen([INTERPRETER, scriptpath] + sys.argv[1:]).wait()
    exit()

# paths
DATA_PATH = '/cs/bd/SequencingData/TranSEQ/'
LOG_PATH = '/cs/bd/TranSEQ/logs'
LOCAL_BIN = '/cs/bd/tools/'
URL_BASE = 'http://www.cs.huji.ac.il/~'

#input
DELIM = ','

# output
UNALIGNED_NAME = 'unaligned'
ANALYSIS_OUTPUT = 'analysis_output'
ANNOTDATA_NAME = 'arrays'
BIGWIG_NAME = 'BIGWIG'
FILTERED_NAME = 'filtered'
TMP_NAME = '.tmp'
NO_BC_NAME = 'no-barcode'
BC_COUNTS_FNAME = 'bc-counts-%s'
SAM_HDR_SUFF = '.sam.hdr'
TMP_BAM_SUFF = '.tmp.bam'
BT_STATS_SUFF = '.bowtie.stats'

# executables
BOWTIE_EXEC = 'bowtie2'
SAMTOOLS_EXEC = 'samtools'
BEDTOOLS_EXEC = 'bedtools'
BG2W_EXEC = LOCAL_BIN+'bedGraphToBigWig'
GC_EXEC = 'genomeCoverageBed'

#error handling
RETRY_INTERVAL = 5  # sec
RETRIALS = 3

#meta
SCER_GENOME_LENGTH_PATH = '/Users/user/Dropbox/workspace/transeq_pipeline/data/sacCer3.sizes'
TTS_MAP = '/cs/wetlab/genomics/scer/annotations/weiner2015_tts.tsv'

STRANDS = {'w': '+', 'c': '-'}


