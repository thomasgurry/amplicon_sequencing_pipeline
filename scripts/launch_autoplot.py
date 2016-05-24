"""

OVERVIEW: 

Script that launches the relevant plotting routines within a given working directory.

"""

from optparse import OptionParser
from AutoPlot import *
from CommLink import *
from SummaryParser import *
    
# Read in arguments for the script                                              

usage = "%prog -i DATASET_DIRECTORY"
parser = OptionParser(usage)
parser.add_option("-i", type="string", dest="dataset_dir")
(options, args) = parser.parse_args()
if( not options.dataset_dir ):
    parser.error("No dataset directory specified.")

# Pipe stdout and stderr to log files
sys.stdout = open(os.path.join(os.getenv("HOME"), 'logs/launch_autoplot_stdout.log'), 'w')
sys.stderr = open(os.path.join(os.getenv("HOME"), 'logs/launch_autoplot_stderr.log'), 'w')

# Get info from summary file
dataset_dir = options.dataset_dir
summary_obj = SummaryParser(os.path.join(dataset_dir, 'summary_file.txt'))
summary_obj.ReadSummaryFile()
dataset_ID = summary_obj.datasetID

# Useful function - returns the intersection of two lists
def intersect(a, b):
    c = [val for val in a if val in b]
    return c

autoplot = AutoPlot(dataset_ID, dataset_dir)


# Plot PiCRUST stuff
working_dir = os.path.join(dataset_dir,'picrust_results') 
for filename in os.listdir(working_dir):
    # Look for L1, L2 and L3 KEGG module files
    split_filename = filename.split('.')
    if(len(intersect(split_filename, ["L1"])) > 0): 
        autoplot.plotKEGGabundances(os.path.join(working_dir,filename), 1)
    elif(len(intersect(split_filename,["L2"])) > 0):
        autoplot.plotKEGGabundances(os.path.join(working_dir,filename), 2)
    elif(len(intersect(split_filename, ["L3"])) > 0):
        autoplot.plotKEGGabundances(os.path.join(working_dir,filename), 3)
