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
            dbotu_flag = 'True'
        try:
            dist = summary_obj.attribute_value_16S['DISTANCE_CRITERIA']
        except:
            dist = 0.1
        try:
            abund = summary_obj.attribute_value_16S['ABUNDANCE_CRITERIA']
        except:
            abund = 10
        try:
            pval = summary_obj.attribute_value_16S['DBOTU_PVAL']
        except:
            pval = 0.0005

    elif amplicon_type == "ITS":
        try:
            dbotu_flag = summary_obj.attribute_value_ITS['DBOTU']
        except:
            dbotu_flag = 'True'
        try:
            dist = summary_obj.attribute_value_ITS['DISTANCE_CRITERIA']
        except:
            dist = 0.1
        try:
            abund = summary_obj.attribute_value_ITS['ABUNDANCE_CRITERIA']
        except:
            abund = 10
        try:
            pval = summary_obj.attribute_value_ITS['DBOTU_PVAL']
        except:
            pval = 0.0005    

    else:
        raise NameError("Incorrect amplicon type specified for dbOTU summary file parsing")

    return dbotu_flag, dist, abund, pval
