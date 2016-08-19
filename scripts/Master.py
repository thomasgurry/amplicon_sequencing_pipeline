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
from Features import *
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

'''
# Load features from ITS processing
feature_dictionary = Features(summary_filename)
print "[[ Feature extraction ]] Loading OTU abundance features for each sample."
feature_dictionary.LoadOTUtable()
print "[[ Feature extraction ]] Loading predicted KEGG module abundances for each sample."
feature_dictionary.LoadPredictedMetagenome()
print "[[ Feature extraction ]] Load L-R phylogenetic features for each sample."
feature_dictionary.LoadPhylogeneticFeatures()

# Pickle feature dictionary
pickled_feature_file = processed_dir + '/' + summary_obj.datasetID + '.features.pk1'
with open(pickled_feature_file, 'wb') as fid:
    pickle.dump(feature_dictionary, fid)

# Navigate to results directory
results_dir = os.path.join('/home/ubuntu/processing_results/', summary_obj.datasetID + '_results')
os.chdir(results_dir)
    
# Set up all files required for User Interface
UI_dir = os.path.join(results_dir, 'UI')
os.system('mkdir ' + UI_dir)
target_directory = os.path.join(results_dir, 'UI/json')
os.system('mkdir ' + target_directory)
UI.write_json_files(summary_obj.attribute_value_16S['OTU_TABLE_CLOSED_REF'], os.path.join(UI_dir, 'json'))
'''
