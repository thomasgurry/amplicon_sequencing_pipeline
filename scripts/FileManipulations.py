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


