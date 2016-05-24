'''

OVERVIEW:  Python module for generating plots for processed data automatically.

'''

import os, sys
import matplotlib.pyplot as plt
from pylab import *

class AutoPlot():
    def __init__(self, datasetID, workingdir):
        # Initialize output directories
        self.datasetID = datasetID
        self.workingdir = os.path.join(workingdir,'plots') 
        self.logfiledir = '/home/ubuntu/logs'
        sys.stdout = open(self.logfiledir + '/listener_stdout.log', 'w')
        sys.stderr = open(self.logfiledir + '/listener_stderr.log', 'w')
        os.system('mkdir ' + self.workingdir)

    def plotKEGGabundances(self, filepath, level):
        # Plots KEGG modules by grouping them to the specified level (levels: 1,2 or 3)
        cmd_lvl = 'summarize_taxa_through_plots.py -i ' + filepath + ' -p /home/ubuntu/qiime_params/qiime_params_L' + str(level) + '.txt -o ' + self.workingdir + '/picrust/plots_at_level' + str(level)
        os.system(cmd_lvl)

    def plotAlphaDiversities(self, alphaDiversityFile, figure_filename):
        # Take an alpha diversity file and create a box plot
        with open(alphaDiversityFile,'r') as fid:
            all_lines = fid.readlines()
            alpha_diversities = [float(line.split()[1]) for line in all_lines[1:]]
            sampleIDs = [line.split()[0] for line in all_lines[1:]]
            figure()
            plt.boxplot(alpha_diversities)
            plt.xlabel('Sample category')
            plt.ylabel('Alpha diversity')
            plt.savefig(figure_filename)



