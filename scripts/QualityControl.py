"""

OVERVIEW: 

Python module for quality control of datasets prior to preprocessing.

"""

import os
import util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Formatting as frmt

def check_length_stats(fastq_in):
    # Returns 25th, 50th, 75th and 95th percentile of read lengths
    iter_seq = util.iter_fsq
    x = []
    counter = 0
    for record in iter_seq(fastq_in):
        sid, seq = record[:2]
        counter = counter + 1
        if(counter > 100000):
            break
        x.append(len(seq))
    x = np.array(x)

    return [np.percentile(x,25), np.percentile(x,50), np.percentile(x,75), np.percentile(x,95)]


def read_length_histogram(raw_sequences_file, path, raw_sequences_filetype='FASTQ'):
    # Creates a histogram of read lengths
    if raw_sequences_filetype == "FASTQ":
        iter_seq = util.iter_fsq
    else:
        iter_seq = util.iter_fst
    x = []
    counter = 0
    for record in iter_seq(raw_sequences_file):
        [sid, seq] = record[:2]
        counter = counter + 1
        if(counter > 100000):
            break
        x.append(len(seq))
    x = np.array(x)
    plt.figure()
    plt.hist(x, 50)
    plt.title('Distribution of amplicon read lengths')
    plt.xlabel('Read length')
    plt.ylabel('Freq')
    plt.savefig(os.path.join(path, 'read_lengths_distribution.png'))


def sample_read_counts(OTU_table, path):
    # Takes as input an OTU table in classic dense format, and the folder to write to, and creates a barchart for the number of sequence reads for each sample.  
    OTU_IDs, sample_IDs, OTU_table = frmt.load_OTU_table_classic(OTU_table)
    readcounts = np.sum(OTU_table, axis=0)
    plt.figure()
    plt.bar(range(len(readcounts)) , readcounts)
    plt.xticks(range(len(readcounts)), sample_IDs, rotation='vertical')
    plt.title('Mean counts per sample = ' + str(int(np.mean(readcounts))))
    plt.xlabel('Sample ID')
    plt.ylabel('Read counts')
    plt.savefig(os.path.join(path, 'sample_read_counts.png'))


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def reads_thrown_out_at_each_step(raw_split_filenames, output_file):
    # Creates a file with % of reads retained after each processing step
    raw_counts = 0
    sb_counts = 0
    pt_counts = 0
    qt_counts = 0
    lt_counts = 0
    fasta_counts = 0

    for filename in raw_split_filenames:
        raw_counts += file_len(filename)
        try:
            sb_counts += file_len(filename + '.sb')
        except:
            pass
        try:
            pt_counts += file_len(filename + '.sb.pt')
        except:
            pass
        try:
            qt_counts += file_len(filename + '.sb.pt.qt')
        except:
            pass
        try:
            lt_counts += file_len(filename + '.sb.pt.qt.lt')
        except:
            pass
        try:
            fasta_counts += file_len(filename + '.sb.pt.qt.lt.fasta')*2
        except:
            pass


        

    with open(output_file, 'a+') as fid:
        line = 'Number of raw reads = ' + str(raw_counts)
        fid.write(line + '\n')
        line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-raw_counts)/float(raw_counts)) + '%'
        fid.write(line + '\n')
        fid.write('\n')

        try:
            line = 'Number of demultiplexed reads = ' + str(sb_counts)
            fid.write(line + '\n')
            line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-sb_counts)/float(raw_counts)) + '%'
            fid.write(line + '\n')
            fid.write('\n')
        except:
            pass
        try:
            line = 'Number of primer-trimmed reads = ' + str(pt_counts)
            fid.write(line + '\n')
            line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-pt_counts)/float(raw_counts)) + '%'
            fid.write(line + '\n')
            fid.write('\n')
        except:
            pass
        try:
            line = 'Number of quality-trimmed reads = ' + str(qt_counts)
            fid.write(line + '\n')
            line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-qt_counts)/float(raw_counts)) + '%'
            fid.write(line + '\n')
            fid.write('\n')
        except:
            pass
        try:
            line = 'Number of length-trimmed reads = ' + str(lt_counts)
            fid.write(line + '\n')
            line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-lt_counts)/float(raw_counts)) + '%'
            fid.write(line + '\n')
            fid.write('\n')
        except:
            pass
        try:
            line = 'Number of FASTA reads left = ' + str(fasta_counts)
            fid.write(line + '\n')
            line = 'Percent of reads left: ' + str(100 - 100*float(raw_counts-fasta_counts)/float(raw_counts)) + '%'
            fid.write(line + '\n')
            fid.write('\n')
        except:
            pass


        
