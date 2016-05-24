"""

OVERVIEW: 

Module to communicate between different servers in the processing layer.

"""

import os, sys
import time
from os.path import basename
from SummaryParser import *

class CommLink():
    def __init__(self, server_type):
        # Takes as input either 'proc', 'picrust' or 'plotting.
        # If the server type is 'picrust', it will launch a listener that will stay on
        # and listen for incoming files.  This should be run from the CRONTAB on a PiCRUST server.
        # If the server type is 'proc', it will launch a listener that will only last until
        # the required files have been received.  This should be run from the 'raw2otu' routine.
        self.procIP = '52.5.178.212'
        self.picrustIP = '52.4.197.174'
        self.plottingIP = '54.152.252.2'
        self.sshkey = '/home/ubuntu/kys/PA_key.pem'
        self.inbox = '/home/ubuntu/inbox'
        self.outbox = '/home/ubuntu/outbox'
        self.logfiledir = '/home/ubuntu/logs'

    def send2picrust(self, filename):
        # Sends filename or folder to PiCRUST server
        if(os.path.isdir(filename)):
            scp_str = 'scp -r -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.picrustIP + ':' + self.inbox + '/.'
        else:
            scp_str = 'scp -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.picrustIP + ':' + self.inbox + '/.'
        os.system(scp_str)
    
    def send2proc(self, filename):
        # Sends filename to processing server
        if(os.path.isdir(filename)):
            scp_str = 'scp -r -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.procIP + ':' + self.inbox + '/.'
        else:
            scp_str = 'scp -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.procIP + ':' + self.inbox + '/.'
        os.system(scp_str)

    def send2plot(self, filename):
        # Sends filename to QIIME/plotting server
        if(os.path.isdir(filename)):
            scp_str = 'scp -r -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.plottingIP + ':' + self.inbox + '/.'
        else:
            scp_str = 'scp -i ' + self.sshkey + ' ' + filename + ' ubuntu@' + self.plottingIP + ':' + self.inbox + '/.'  
        os.system(scp_str)

    def intersect(self, a, b):
        # Returns the intersection of two lists
        c = [val for val in a if val in b]
        return c

    def launch_proc_listener(self, outgoing_directory, expected_results_directory):
        # Launches a listener that sends a folder to PiCRUST server and waits until the expected directory is received.        
        
        # Send all files to PiCRUST server
        self.send2picrust(outgoing_directory)
        summary_obj = SummaryParser(os.path.join(outgoing_directory, 'summary_file.txt'))
        summary_obj.ReadSummaryFile()
        datasetID = summary_obj.datasetID

        # Send file telling proc server the data transfer is complete
        with open(datasetID + '_1_3_transfer.txt', 'w') as fid:
            fid.write('Complete.')
        self.send2picrust(datasetID + '_1_3_transfer.txt')

        # Waits to get results back before proceeding.
        inbox_files = os.listdir(self.inbox)
        transfer_complete = False
        while(len(self.intersect(inbox_files, expected_results_directory)) == 0 and transfer_complete == False):
            # Sleep 2 seconds and check again
            time.sleep(2)
            inbox_files = os.listdir(self.inbox)        
            for directory in inbox_files: 
                # Check summary file in directory, and check for 'transfer complete file'
                try:
                    summary_obj = SummaryParser(os.path.join(self.inbox, directory, 'summary_file.txt'))
                    summary_obj.ReadSummaryFile()
                    datasetID = summary_obj.datasetID
                    transfer_complete = os.path.isfile(os.path.join(self.inbox, datasetID + '_3_3_transfer.txt'))
                    if(transfer_complete == False):
                        continue                        
                    else:
                        break
                except:
                    continue
                if(transfer_complete == True):
                    break
        # Delete 'transfer complete file' and close listener
        os.system('rm ' + os.path.join(self.inbox, datasetID + '_3_3_transfer.txt'))
        return 1

        # Wait for transfer to complete (receives a file with name 'datasetID_transfer.txt')
        #        while(len(self.intersect(inbox_files, transfer_complete_file)) == 0):
            # Sleep 2 seconds and check again
            #            time.sleep(2)
            #            inbox_files = os.listdir(self.inbox)        

    def launch_picrust_listener(self, inbox_folder):
        # Launches a listener for the PiCRUST server, that will listen
        # indefinitely for incoming dataset directories in the specified inbox folder, perform
        # metagenomic prediction using PiCRUST, and then place the relevant 
        # files in the outbox before emptying the inbox.        

        # Pipe stdout and stderr to log files
        sys.stdout = open(self.logfiledir + '/listener_stdout.log', 'w')
        sys.stderr = open(self.logfiledir + '/listener_stderr.log', 'w')

        # Folder where work takes place
        work_folder = '/home/ubuntu/picrust_working'
        
        while(1):
            # Check every 2 seconds.
            time.sleep(2)
            # Get contents of inbox folder specified
            folders_in_inbox = [d for d in os.listdir(inbox_folder) if os.path.isdir(os.path.join(inbox_folder, d))]
            if len(folders_in_inbox) > 0:
                for directory in folders_in_inbox:                    
                    # Check summary file in directory, and check for 'transfer complete file'
                    try:
                        summary_obj = SummaryParser(os.path.join(inbox_folder, directory, 'summary_file.txt'))
                        summary_obj.ReadSummaryFile()
                        datasetID = summary_obj.datasetID
                        transfer_complete = os.path.isfile(os.path.join(inbox_folder, datasetID + '_1_3_transfer.txt'))
                        if(transfer_complete == False):
                            continue
                    except:
                        continue

                    # Make working directory and move file there
                    os.system('mv ' + inbox_folder + '/' + directory + ' ' + work_folder + '/.')
                    os.system('rm ' + os.path.join(inbox_folder, datasetID + '_1_3_transfer.txt'))
                    picrust_working_dir = work_folder + '/' + directory

                    # Parse summary file
                    summary_obj = SummaryParser(picrust_working_dir + '/summary_file.txt')
                    summary_obj.ReadSummaryFile()
                    datasetID = summary_obj.datasetID
                    OTU_table_file_closed_ref = picrust_working_dir + '/' + summary_obj.attribute_value_16S['OTU_TABLE_CLOSED_REF']
                    OTU_table_file_open_ref = picrust_working_dir + '/' + summary_obj.attribute_value_16S['OTU_TABLE_OPEN_REF']

                    # Initialize PiCRUST-related filenames
                    qiime_ready_file_closed = OTU_table_file_closed_ref + '.qiime_ready.txt'
                    norm_filename_closed = OTU_table_file_closed_ref + '.norm.biom'
                    metagenome_predictions = picrust_working_dir + '/' + datasetID + '.metapred.biom'
                    nsti_file = picrust_working_dir + '/' + datasetID + '.nsti_values.tab'
                    level_3_predictions = picrust_working_dir + '/' + datasetID + '.L3.biom'
                    level_2_predictions = picrust_working_dir + '/' + datasetID + '.L2.biom'
                    level_1_predictions = picrust_working_dir + '/' + datasetID + '.L1.biom'
                    otu_table_closed_with_metadata = picrust_working_dir + '/' + datasetID + '.otu_table_closed_with_metadata.biom'
                    otu_table_open_with_metadata = picrust_working_dir + '/' + datasetID + '.otu_table_open_with_metadata.biom'

                    # Convert OTU tables to BIOM format and add metadata for each sample if metadata exists
                    os.system('biom convert -i ' + OTU_table_file_closed_ref + ' -o ' + picrust_working_dir + '/' + datasetID + '.otu_table_closed.biom --table-type="OTU table"')
                    os.system('biom convert -i ' + OTU_table_file_open_ref + ' -o ' + picrust_working_dir + '/' + datasetID + '.otu_table_open.biom --table-type="OTU table"')
                    try:
                        metadata_file = picrust_working_dir + '/' + summary_obj.attribute_value_16S['METADATA_FILE']
                        os.system('biom add-metadata -i ' + picrust_working_dir + '/' + datasetID + '.otu_table_closed.biom -m ' + metadata_file + ' -o ' + otu_table_closed_with_metadata)
                        os.system('biom add-metadata -i ' + picrust_working_dir + '/' + datasetID + '.otu_table_open.biom -m ' + metadata_file + ' -o ' + otu_table_open_with_metadata)
                    except:
                        metadata_file = None
                        print "No metadata file specified in the summary file.  Proceeding without metadata."
                        otu_table_closed_with_metadata = os.path.join(picrust_working_dir, datasetID + '.otu_table_closed.biom')

                    os.system('biom convert -b -i ' + otu_table_closed_with_metadata + ' -o ' + qiime_ready_file_closed + ' --table-type="OTU table"')

                    # PiCRUST pipeline
                    # Step 1: normalize OTU counts by copy number
                    os.system('normalize_by_copy_number.py -f -i ' + qiime_ready_file_closed + ' -o ' + norm_filename_closed)

                    # Step 2: predict metagenomes
                    os.system('predict_metagenomes.py -i ' + norm_filename_closed + ' -o ' + metagenome_predictions + ' -a ' + nsti_file)

                    # Step 3: collapse predictions into higher categories e.g. KOs into KEGG pathways
                    cmd_str3 = 'categorize_by_function.py -i ' + metagenome_predictions + ' -c KEGG_Pathways -l 3 -o ' + level_3_predictions
                    cmd_str2 = 'categorize_by_function.py -i ' + metagenome_predictions + ' -c KEGG_Pathways -l 2 -o ' + level_2_predictions
                    cmd_str1 = 'categorize_by_function.py -i ' + metagenome_predictions + ' -c KEGG_Pathways -l 1 -o ' + level_1_predictions
                    os.system(cmd_str3)
                    os.system(cmd_str2)
                    os.system(cmd_str1)

                    # Step 4: move results to a single folder and send to plotting server to parse results
                    picrust_results_folder = picrust_working_dir + '/' + 'picrust_results'
                    os.system('mkdir ' + picrust_results_folder)
                    os.system('mv ' + metagenome_predictions + ' ' + picrust_results_folder + '/.')
                    summary_obj.attribute_value_16S['METAGENOME_PREDICTIONS'] = os.path.join('picrust_results', basename(metagenome_predictions))
                    # move KEGG predictions
                    os.system('mv ' + level_3_predictions + ' ' + picrust_results_folder + '/.')
                    os.system('mv ' + level_2_predictions + ' ' + picrust_results_folder + '/.')
                    os.system('mv ' + level_1_predictions + ' ' + picrust_results_folder + '/.')
                    summary_obj.attribute_value_16S['KEGG_L1'] = os.path.join('picrust_results', basename(level_1_predictions))
                    summary_obj.attribute_value_16S['KEGG_L2'] = os.path.join('picrust_results', basename(level_2_predictions))
                    summary_obj.attribute_value_16S['KEGG_L3'] = os.path.join('picrust_results', basename(level_3_predictions))
                    summary_obj.WriteSummaryFile()
                    self.send2plot(picrust_working_dir)

                    # Send file telling proc server the data transfer is complete
                    with open(datasetID + '_2_3_transfer.txt', 'w') as fid:
                        fid.write('Complete.')
                    self.send2plot(datasetID + '_2_3_transfer.txt')


    def launch_plotting_listener(self, inbox_folder):
        # Launches a listener for the QIIME/plotting server.

        # Pipe stdout and stderr to log files
        sys.stdout = open(self.logfiledir + '/listener_stdout.log', 'w')
        sys.stderr = open(self.logfiledir + '/listener_stderr.log', 'w')

        # Folder where work takes place
        work_folder = '/home/ubuntu/plot_working'
        
        while(1):
            # Check every 2 seconds.
            time.sleep(2)
            # Get folders in inbox folder specified
            folders_in_inbox = [d for d in os.listdir(inbox_folder) if os.path.isdir(os.path.join(inbox_folder, d))]
            if len(folders_in_inbox) > 0:
                for folder in folders_in_inbox:                    
                    # Check summary file in directory, and check for 'transfer complete file'
                    try:
                        summary_obj = SummaryParser(os.path.join(inbox_folder, folder, 'summary_file.txt'))
                        summary_obj.ReadSummaryFile()
                        datasetID = summary_obj.datasetID
                        transfer_complete = os.path.isfile(os.path.join(inbox_folder, datasetID + '_2_3_transfer.txt'))
                        if(transfer_complete == False):
                            continue
                    except:
                        continue

                    # Move the folder to the working directory
                    os.system('mv ' + os.path.join(inbox_folder, folder) + ' ' + work_folder + '/.')
                    os.system('rm ' + os.path.join(inbox_folder, datasetID + '_2_3_transfer.txt'))
                    working_dir = os.path.join(work_folder, folder)

                    # Parse summary file
                    summary_obj = SummaryParser(os.path.join(working_dir, 'summary_file.txt'))
                    summary_obj.ReadSummaryFile()
                    datasetID = summary_obj.datasetID
                    OTU_sequences_open_ref = os.path.join(working_dir, summary_obj.attribute_value_16S['OTU_SEQUENCES_FASTA'])
                    OTU_table_file_closed_ref = os.path.join(working_dir, summary_obj.attribute_value_16S['OTU_TABLE_CLOSED_REF'])
                    OTU_table_file_open_ref = os.path.join(working_dir, summary_obj.attribute_value_16S['OTU_TABLE_OPEN_REF'])
                    metadata_file = os.path.join(working_dir, summary_obj.attribute_value_16S['METADATA_FILE'])

                    # Compute phylogenetic tree for open reference OTUs
                    phylogenetic_tree_file = os.path.join(working_dir, 'phylogenetic_tree_open_reference.tre')
                    os.system('make_phylogeny.py -i ' + OTU_sequences_open_ref + ' -o ' + phylogenetic_tree_file)
                    summary_obj.attribute_value_16S['PHYLOGENETIC_TREE_OPEN_REF'] = basename(phylogenetic_tree_file)

                    # Compute alpha rarefaction curves until a sequencing depth equal to the minimum number of reads found in the samples
                    with open(OTU_table_file_open_ref, 'r') as fid:
                        otu_table_data = fid.readlines()                                      
                        firstrow = otu_table_data[0].split('\t')                                
                        OTU_labels = firstrow[1:]                                               
                        OTU_labels[len(OTU_labels)-1] = OTU_labels[len(OTU_labels)-1].rstrip('\n')
                        sample_labels = [otu_table_data[i].split('\t')[0] for i in range(1,len(otu_table_data))]                                                               
                        nOTUs = len(OTU_labels)                                                 
                        nSamples = len(sample_labels)                                           
                        # Load OTU table row major order                                        
                        OTU_table = np.zeros((nSamples, nOTUs))                             
                        for i in range(1,nSamples+1):                                           
                            OTU_table[i-1,:] = otu_table_data[i].split('\t')[1:]            
                        min_read_count = int(np.min(np.sum(OTU_table, axis=0)))
                        
                    alpha_rarefaction_dir = os.path.join(working_dir, 'alpha_rarefaction')
                    os.system('alpha_rarefaction.py -i ' + OTU_table_file_open_ref + ' -m ' + metadata_file + ' -o ' + alpha_rarefaction_dir + ' -t ' + phylogenetic_tree_file + ' -e ' + str(min_read_count))

                    # Compare alpha diversities
                    split_by = 'DiseaseState'
                    alpha_div_comparisons_directory = os.path.join(working_dir, 'plots', 'alpha_diversity')
                    alpha_div_file = os.path.join(alpha_rarefaction_dir,'alpha_div_collated','chao1.txt')
                    os.system('compare_alpha_diversity.py -i ' + alpha_div_file + ' -m ' + metadata_file + ' -c ' + split_by + ' -o ' + alpha_div_comparisons_directory)

                    # Launch AutoPlot routines
                    os.system('python ~/scripts/launch_autoplot.py -i ' + working_dir)

                    # Update summary file
                    summary_obj.WriteSummaryFile()

                    # Send back to proc server
                    self.send2proc(working_dir)
                    
                    # Send file telling proc server the data transfer is complete
                    with open(datasetID + '_3_3_transfer.txt', 'w') as fid:
                        fid.write('Complete.')
                    self.send2proc(datasetID + '_3_3_transfer.txt')
