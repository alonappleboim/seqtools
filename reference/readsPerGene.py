# Input - directory with BAM files ,interval type (ORF,TSS-TTS),  output path (directory/filename.csv)
# interval types options: ORF, TTS-250, TTS-350, TTS-500, TTS-750, TTS-1000...
#   (as for now it is half-parsed as a string)
# Output - csv with counts (not normalized), and report for the counted reads.

import os
import sys
import csv
import pysam

TSS_TTS_path = "/cs/bd/SequencingData/RNA-MAAYAN-ADI/scripts/SGD_TSS.csv"


def handlebam(fname, type):
    """
    For a single bam file, it parses the sgd_tss positions file, and for each gene in this file
    it looks for reads that are fetched to this position.
    """
    print 'Current file: %s' % fname
    bamfile = pysam.AlignmentFile(fname, "rb")
    totalreads = (bamfile.count())

    count = []
    # parse the genes positions file
    with open(TSS_TTS_path) as input_genes:
        for lines in input_genes:
            lines = lines.split('\r')
        for line in lines[1:]:
            line = line.split('\t')
            # chr17 is not indexed in our bam file and throws an error, so do not try to fetch it:
            if line[2] == 'XVII' or line[2] == 'chrXVII':
                continue

            intervals = getIntervals(type, line)
            s = min(int(intervals[0]), int(intervals[1]))
            e = max(int(intervals[0]), int(intervals[1]))
            rev = int(intervals[0]) > int(intervals[1])

            #format the chr name IMPORTANT - this change is needed only in part of the bam file:
            #if not line[2].startswith('chr'):
                #line[2] = 'chr' + line[2]
    
                        
            # look for reads in this position (interval). If there are no reads, continue
            try: 
                tmp = bamfile.fetch(line[2], s, e)
            except ValueError:
                print "ValueError"
                print line[2]
                print s 
                print e
                print "========================"
                continue
            
            # count only reads of the same strand as the gene
            cnt = 0
            for r in tmp:
                if rev == r.is_reverse:
                    cnt += 1

            count.append(float(cnt))

    counted_reads = sum(count)

    return count, totalreads, counted_reads


def getIntervals(type, line):
    intervals = [0, 0]
    rev = int(line[3]) > int(line[4])
    if type == 'ORF':
        intervals[0] = line[3]
        intervals[1] = line[4]
    else:
        # Check if TTS/TSS fit the ORF
        if line[5] != "NaN":
            # TSS after start:
            if not rev and int(line[3]) < int(line[5]):
                line[5] = "NaN"
            if rev and int(line[3]) > int(line[5]):
                line[5] = "NaN"
        if line[6] != "NaN":
            # TTS before end:
            if not rev and int(line[4]) > int(line[6]):
                line[6] = "NaN"
            if rev and int(line[4]) < int(line[6]):
                line[6] = "NaN"

        # If TSS/TTS don't fit the direction of the ORF
        if line[5] != "NaN" and line[6] != "NaN":
            rev2 = int(line[5]) > int(line[6])
            if rev != rev2:
                line[5] = "NaN"
                line[6] = "NaN"
        # fill the start of the interval with TSS if possible, else with 'start ORF':
        if line[5] == "NaN":
            intervals[0] = int(line[3])
        else:
            intervals[0] = int(line[5])

        # If there is no TTS data we take end of ORF+300 otherwise TTS+200
        if line[6] == "NaN" or int(line[6]) == int(line[4]):
            gene_end = int(line[4])
            if rev:
                intervals[1] = gene_end-300
            else:
                intervals[1] = gene_end+300
        else:
            gene_end = int(line[6])
            if rev:
                intervals[1] = gene_end-200
            else:
                intervals[1] = gene_end+200

        # update the start position according to input:
        factor = int(type.split('-')[1])
        intervals = updateIntervals(intervals, rev, factor, gene_end)

    if intervals[0] < 0:
        intervals[0] = 0
    if intervals[1] < 0:
        intervals[1] = 0

    return intervals


def updateIntervals(intervals, rev, factor, gene_end):
    # update the interval to start "factor" before the gene_end, if it doesn't pass the start.
    if not rev:
        tmp = gene_end-factor
        if tmp > int(intervals[0]):
            intervals[0] = tmp
    else:
        tmp = gene_end+factor
        if tmp < int(intervals[0]):
            intervals[0] = tmp

    return intervals


def main(args):
    print "Type: %s" % args[2]
    dirname = os.path.abspath(args[1])
    filenames = []
    # take all bam files from the input dir:
    for fname in os.listdir(dirname):
        if not(fname.startswith('.')) and fname.endswith('.bam'):
            filenames.append('/'.join([dirname, fname]))

    with open(TSS_TTS_path) as input_genes:
        for lines in input_genes:
            lines = lines.split('\r')
        genesnum = len(lines)-1

    # initialize matrix for genes X samples
    datamat = [[0 for x in range(len(filenames)+6)] for x in range(genesnum+1)]

    datamat[0][0] = "ORF"
    datamat[0][1] = "GeneName"
    datamat[0][2] = "Chromosome"
    datamat[0][3] = "Start"
    datamat[0][4] = "End"
    datamat[0][5] = "Reverse"

    # insert genes details to the first 6 cols of datamat
    ind = 1
    for line in lines[1:]:
        line = line.split('\t')
        line[0] = line[0].replace("\"", "")
        datamat[ind][0] = line[0]
        if not line[1]:
            datamat[ind][1] = line[0]
        else:
            datamat[ind][1] = line[1].replace(',', ';')

        datamat[ind][2] = line[2]
        intervals = getIntervals(args[2], line)
        datamat[ind][3] = int(intervals[0])
        datamat[ind][4] = int(intervals[1])
        datamat[ind][5] = datamat[ind][3] > datamat[ind][4]
        ind += 1

    report = [['sampleName', 'totalReads', 'countedReads']]

    # fill table with read counts per each sample - NOT NORMALIZED:
    ind = 0
    for fname in filenames:

        # if ind > 16:
        #     continue
        # sample_name = 'xxx'

        # insert the col name for each sample: ###IMPORTANT - change for each experiment###
        tmp_name = os.path.basename(fname).split('_')
        if tmp_name[3] == 'A' or tmp_name[3] == 'B':
            sample_name = '_'.join(tmp_name[:4])
        else:
            sample_name = '_'.join(tmp_name[:3])

        datamat[0][ind+6] = sample_name
        ####datamat[0][ind+6] = os.path.basename(fname).split('.')[0]

        # count reads per genes for this sample and insert to the table:
        out, total_reads, counted_reads = handlebam(fname, args[2])
        for i in range(len(out)-1):
            datamat[i+1][ind+6] = (out[i])
        report.append([sample_name, total_reads, counted_reads])

        ind += 1

    # sort the cols by alphabet order: for now it is done in the normalization script :
    # normalize_order_readcounts.py, that runs on this script's output

    # save output file
    outfilename = os.path.abspath(args[3])
    outfile = open(outfilename, 'w')
    csvf = csv.writer(outfile)
    csvf.writerows(datamat)
    outfile.close()

    # save the read counts report
    reportname = '.'.join(outfilename.split('.')[:-1]) + '_report.csv'
    reportfile = open(reportname, 'w')
    csvf2 = csv.writer(reportfile)
    csvf2.writerows(report)
    reportfile.close()

if __name__ == '__main__':
    main(sys.argv)