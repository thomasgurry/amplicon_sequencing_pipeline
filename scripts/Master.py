"""

OVERVIEW: 

Master Script that checks a given dataset directory for any processing requests.

Currently only 16S processing supported.

"""

import numpy as np
import os
import os.path
import sys
import pickle
from optparse import OptionParser
from SummaryParser import *
import UserInterface as UI
import QualityControl as QC

# Read in arguments
usage = "%prog -i DATASET_DIR"
parser = OptionParser(usage)
parser.add_option("-i", "--datadir", type="string", dest="datadir")
(options, args) = parser.parse_args()

if( not options.datadir ):
    parser.error("No directory specified for the data.")


# Pipe stdout and stderr to logfiles in the new directory
working_directory = options.datadir
sys.stdout = open(os.path.join(working_directory, 'stdout_master.log'),'w')
sys.stderr = open(os.path.join(working_directory, 'stderr_master.log'),'w')

# Parse summary file
summary_filename = os.path.join(options.datadir, 'summary_file.txt')
summary_obj = SummaryParser(summary_filename)
summary_obj.ReadSummaryFile()

# Check if 16S processing is required, and if so, launch the processing pipeline.
if(summary_obj.attribute_value_16S['PROCESSED'] == 'True'):
    print "[[ 16S processing ]] Processing already complete."
elif(summary_obj.attribute_value_16S['PROCESSED'] == 'False'):

    flags = ''

    # Check if multiple FASTQ files are provided.  If so, assuming they are demultiplexed (each file is a particular sample/barcode).
    try:
        raw_data_summary_file = summary_obj.attribute_value_16S['RAW_FASTQ_FILES']
        flags = flags + ' -m True '
    except:
        try:
            raw_data_summary_file = summary_obj.attribute_value_16S['RAW_FASTA_FILES']
            flags = flags + ' -m True '
        except:
            flags = flags + ' -m False '
    # Check if paired end reads need to be merged
    try:
        paired_ends = summary_obj.attribute_value_16S['MERGE_PAIRS']
        flags = flags + ' -r ' + paired_ends
    except:
        pass

    # Check if output directory is specified
    try:
        outdir = summary_obj.attribute_value_16S['OUTDIR'] 
        flags = flags + ' -o ' + outdir
    except:
        pass
        
    # Check if primers have been removed
    if summary_obj.attribute_value_16S['PRIMERS_FILE'] == 'None':
        flags = flags + ' -p True'
        print "[[ 16S processing ]] No primers file.  Assuming primers have been trimmed."
    if summary_obj.attribute_value_16S['BARCODES_MAP'] == 'None':
        flags = flags + ' -b True'
        print "[[ 16S processing ]] No barcodes map.  Assuming sequences have been demultiplexed and relabeled with sample IDs."
    print "[[ 16S processing ]] Processing required.  Generating OTU tables."
    raw2otu_cmd = 'python ~/scripts/raw2otu.py -i ' + working_directory + flags 
    os.system(raw2otu_cmd)

    # Check summary file again
    summary_obj.ReadSummaryFile()
    if(summary_obj.attribute_value_16S['PROCESSED'] == 'False'):
        print "[[ 16S processing ]] ERROR: Failed to process 16S data"
        raise NameError('ERROR: Failed to process 16S data.')    
    elif(summary_ob.attribute_value_16S['PROCESSED'] == 'True'):        
        print "[[ 16S processing ]] Processing complete."
    else:
        print "[[ 16S processing ]] Unknown case.  Proceeding..."

# Check if ITS processing is required, and if so, launch the processing pipeline.
if(summary_obj.attribute_value_ITS['PROCESSED'] == 'True'):
    print "[[ ITS processing ]] Processing already complete."
elif(summary_obj.attribute_value_ITS['PROCESSED'] == 'False'):

    flags = ''

    # Check if multiple FASTQ files are provided.  If so, assuming they are demultiplexed (each file is a particular sample/barcode).
    try:
        raw_data_summary_file = summary_obj.attribute_value_ITS['RAW_FASTQ_FILES']
        flags = flags + ' -m True '
    except:
        try:
            raw_data_summary_file = summary_obj.attribute_value_ITS['RAW_FASTA_FILES']
            flags = flags + ' -m True '
        except:
            flags = flags + ' -m False '

    # Check if primers have been removed
    if summary_obj.attribute_value_ITS['PRIMERS_FILE'] == 'None':
        flags = flags + ' -p True'
        print "[[ ITS processing ]] No primers file.  Assuming primers have been trimmed."
    if summary_obj.attribute_value_ITS['BARCODES_MAP'] == 'None':
        flags = flags + ' -b True'
        print "[[ ITS processing ]] No barcodes map.  Assuming sequences have been demultiplexed and relabeled with sample IDs."
    print "[[ ITS processing ]] Processing required.  Generating OTU tables."
    raw2otu_cmd = 'python ~/scripts/raw2otu.py -i ' + working_directory + flags 
    os.system(raw2otu_cmd)

    # Check summary file again
    summary_obj.ReadSummaryFile()
    if(summary_obj.attribute_value_ITS['PROCESSED'] == 'False'):
        print "[[ ITS processing ]] ERROR: Failed to process ITS data"
        raise NameError('ERROR: Failed to process ITS data.')    
    elif(summary_ob.attribute_value_ITS['PROCESSED'] == 'True'):        
        print "[[ ITS processing ]] Processing complete."
    else:
        print "[[ ITS processing ]] Unknown case.  Proceeding..."

else:
    print "[[ ITS processing ]] No processing request specified."

