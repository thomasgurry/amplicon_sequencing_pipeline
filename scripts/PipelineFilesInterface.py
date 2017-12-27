"""

Helper functions to move files around.

"""

import os
from string import ascii_lowercase

def split_file_and_return_names(raw_data_file):
    """
    Splits an input raw data file into files of different
    line lengths. Returns the created files.
    """

    rawfilesize = os.path.getsize(raw_data_file)

    # Step 1.1 - split file into 1000000 line (~100Mb) chunks
    # Note: the 'split' command takes the raw_data_file, splits it into separate files
    # and puts those files in the current directory (i.e. working_directory)
    if(rawfilesize < 2e8):
        os.system('split -l 100000 ' + raw_data_file)
    else:
        os.system('split -l 1000000 ' + raw_data_file)

    # Step 1.2 - get split filenames
    split_filenames = []
    for c1 in ascii_lowercase:
        for c2 in ascii_lowercase:
            filename = 'x'+c1+c2
            if(os.path.isfile(filename)):
                split_filenames.append(filename)
    '''
    # Claire's note: I'm pretty sure split -l 10000 raw_data_file will always make
    # at least xaa, and so len(split_filenames) should never be zero...
    if len(split_filenames) == 0:
        split_filenames = [raw_data_file]
    raw_filenames = split_filenames
    '''
    return split_filenames


##########################################################################
##% Functions to parse specific summary_file attributes for 
##% different raw2otu.py processes
##########################################################################
def parse_input_files(options, summary_obj, amplicon_type):
    """
    Parses the raw data file input and returns the file name
    and file type.

    Parameters
    ----------
    options                 The options that were passed to raw2otu.py.
                            Should have options.input_dir
    summary_obj             SummaryParser object
    amplicon_type           '16S' or 'ITS'
    
    Returns
    -------
    raw_data_file           raw data file path (i.e. input_dir/raw_file)
    raw_data_summary_file   file path to file with raw data to sample ID map
    raw_file_type           'FASTQ' or 'FASTA'
    barcodes_map            barcodes map file
    primers_files           primers file
    """

    raw_data_file = None
    raw_data_summary_file = None

    # Extract file locations
    if amplicon_type == '16S':
        primers_file = os.path.join(options.input_dir, summary_obj.attribute_value_16S['PRIMERS_FILE'])
        barcodes_map = os.path.join(options.input_dir, summary_obj.attribute_value_16S['BARCODES_MAP'])
        try:
            raw_data_file = os.path.join(options.input_dir, summary_obj.attribute_value_16S['RAW_FASTQ_FILE'])
            raw_file_type = 'FASTQ'
        except:
            print("No single raw FASTQ file found.  Checking for raw FASTA.")
            try:
                raw_data_file = os.path.join(options.input_dir, summary_obj.attribute_value_16S['RAW_FASTA_FILE'])
                raw_file_type = 'FASTA'
            except:
                print("No single raw FASTA file found either.  Checking for multiple files.")
                try:
                    raw_data_summary_file = os.path.join(options.input_dir, summary_obj.attribute_value_16S['RAW_FASTQ_FILES'])
                    raw_file_type = 'FASTQ'
                except:
                    print("No filename of multiple raw FASTQs map provided.  Check contents of your raw data and summary file.")
                    try:
                        raw_data_summary_file = os.path.join(options.input_dir, summary_obj.attribute_value_16S['RAW_FASTA_FILES'])
                        raw_file_type = 'FASTA'
                    except:
                        print("No filename of multiple raw FASTAs map provided.  Check contents of your raw data and summary file.")
                        raise NameError("Unable to retrieve raw sequencing files.")

    elif amplicon_type == 'ITS':
        primers_file = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['PRIMERS_FILE'])
        barcodes_map = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['BARCODES_MAP'])
        try:
            raw_data_file = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['RAW_FASTQ_FILE'])
            raw_file_type = 'FASTQ'
        except:
            print("No single raw FASTQ file found.  Checking for raw FASTA.")
            try:
                raw_data_file = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['RAW_FASTA_FILE'])
                raw_file_type = 'FASTA'
            except:
                print("No single raw FASTA file found either.  Checking for multiple files.")
                try:
                    raw_data_summary_file = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['RAW_FASTQ_FILES'])
                    raw_file_type = 'FASTQ'
                except:
                    print("No filename of multiple raw FASTQs map provided.  Check contents of your raw data and summary file.")
                    raise NameError("Unable to retrieve raw sequencing files.")
                    try:
                        raw_data_summary_file = os.path.join(options.input_dir, summary_obj.attribute_value_ITS['RAW_FASTA_FILES'])
                        raw_file_type = 'FASTA'
                    except:
                        print("No filename of multiple raw FASTAs map provided.  Check contents of your raw data and summary file.")
                        raise NameError("Unable to retrieve raw sequencing files.")

    return primers_file, barcodes_map, raw_data_file, raw_data_summary_file, raw_file_type

def parse_barcodes_parameters(summary_obj, amplicon_type):
    """
    Parses parameters used in splitting by barcodes.
    
    Parameters
    ----------
    summary_obj         SummaryParser object
    amplicon_type       '16S' or 'ITS'

    Returns
    -------
    mode                 str,  either '1', '2', or '3' - indicating
                         barcodes mode to pass to 
                         2.split_by_barcodes.py
    index_file           str, If mode '3' is given, the INDEX_FILE 
                         location, otherwise returns None.
    index_file_format    str (default = 'fastq'). 'fastq', 'fasta',
                         or 'tab' - for use in 2.split_by_barcodes.py
    """
    if amplicon_type == '16S':
        if 'BARCODES_MODE' in summary_obj.attribute_value_16S:
            mode = summary_obj.attribute_value_16S['BARCODES_MODE']
        else:
            mode = '2'
        # If barcodes are in index file, get those parameters
        if mode == '3':
            try:
                index_file = summary_obj.attribute_value_16S['INDEX_FILE']
            except:
                raise NameError("Barcodes mode 3 specified (barcodes are in index file), but no INDEX_FILE provided.")
            try:
                index_file_format = summary_obj.attribute_value_16S['INDEX_FORMAT']
            except:
                index_file_format = 'fastq'
        else:
            index_file = None
            index_file_format = None
    elif amplicon_type == 'ITS':
        if 'BARCODES_MODE' in summary_obj.attribute_value_ITS:
            mode = summary_obj.attribute_value_ITS['BARCODES_MODE']
        else:
            mode = '2'
        # If barcodes are in index file, get those parameters
        if mode == '3':
            try:
                index_file = summary_obj.attribute_value_ITS['INDEX_FILE']
            except:
                raise NameError("Barcodes mode 3 specified (barcodes are in index file), but no INDEX_FILE provided.")
            try:
                index_file_format = summary_obj.attribute_value_ITS['INDEX_FORMAT']
            except:
                index_file_format = 'fastq'
        else:
            index_file = None
            index_file_format = None
    return mode, index_file, index_file_format

def parse_dbotu_parameters(summary_obj, amplicon_type):
    """
    Parses summary file for dbOTU options.

    Parameters
    ----------
    summary_obj         SummaryParser object
    amplicon_type       '16S' or 'ITS'

    Returns
    -------
    dist                max sequence dissimilarity (default = 0.1)
    abund               minimum fold difference for comparing two OTUs (0=no abundance criterion; default 10.0)
    pval                minimum p-value for merging OTUs (default: 0.0005)
    """

    if amplicon_type == "16S":
        try:
            dbotu_flag = summary_obj.attribute_value_16S['DBOTU']
        except:
            dbotu_flag = 'False'
        try:
            dist = summary_obj.attribute_value_16S['DISTANCE_CRITERIA']
        except:
            dist = 0.1
        try:
            abund = summary_obj.attribute_value_16S['ABUNDANCE_CRITERIA']
        except:
            abund = 10.0
        try:
            pval = summary_obj.attribute_value_16S['DBOTU_PVAL']
        except:
            pval = 0.0005

    elif amplicon_type == "ITS":
        try:
            dbotu_flag = summary_obj.attribute_value_ITS['DBOTU']
        except:
            dbotu_flag = 'False'
        try:
            dist = summary_obj.attribute_value_ITS['DISTANCE_CRITERIA']
        except:
            dist = 0.1
        try:
            abund = summary_obj.attribute_value_ITS['ABUNDANCE_CRITERIA']
        except:
            abund = 10.0
        try:
            pval = summary_obj.attribute_value_ITS['DBOTU_PVAL']
        except:
            pval = 0.0005    

    else:
        raise NameError("Incorrect amplicon type specified for dbOTU summary file parsing")

    return dbotu_flag, dist, abund, pval
