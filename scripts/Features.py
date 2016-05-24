"""

OVERVIEW: 

Script to convert processing results into features for use in the machine learning/analytics layer.

"""

import os, sys
from SummaryParser import *
import Formatting as frmt
import Phylogeny as phylo

class Features():
    def __init__(self, summary_file=0):    
        if (summary_file != 0):
            # Parse summary file
            self.summary_file = summary_file
            summary_obj = SummaryParser(summary_file)
            summary_obj.ReadSummaryFile()
            self.dataset_ID = summary_obj.datasetID
            try:
                self.otu_table_closed = summary_obj.attribute_value_16S['OTU_TABLE_CLOSED_REF']
                self.otu_table_open = summary_obj.attribute_value_16S['OTU_TABLE_OPEN_REF']
                # Load sample IDs from OTU table -
                # NOTE: CHANGE THIS TO LOAD SAMPLE IDs FROM STANDARD METADATA FILE
                with open(self.otu_table_closed, 'r') as fid:
                    otu_table_data = fid.readlines()
                    firstrow = otu_table_data[0].split('\t')
                    sample_labels = firstrow[1:]
                    sample_labels[len(sample_labels)-1] = sample_labels[len(sample_labels)-1].rstrip('\n')
            
                # Initialize feature vectors and feature descriptions
                self.feature_descriptions = {}
                self.feature_vectors = {}
                for sample_ID in sample_labels:
                    self.feature_vectors[sample_ID] = []

            except:
                print "Error - failure to load OTU table filename from summary file.  Check summary file contents."
            try:
                self.metapred_table = summary_obj.attribute_value_16S['METAGENOME_PREDICTIONS']
            except:
                print "Error - failure to load metagenome predictions from summary file.  Check summary file contents."
            try:
                self.phylogenetic_tree_open = summary_obj.attribute_value_16S['PHYLOGENETIC_TREE_OPEN_REF']
            except:
                print "Error - failure to load phylogenetic tree from summary file.  Check summary file contents."

    def LoadOTUtable(self, OTU_table_filename=0):
        # Loads an OTU table into the feature vectors for each sample.  Assumes the  OTU table is in the final form you want, i.e. normalized as desired, etc.  Will default to the closed reference OTU table in the summary file unless otherwise specified.

        if OTU_table_filename == 0:
            OTU_table_filename = self.otu_table_closed
        with open(OTU_table_filename,'r') as fidin:
            otu_table_data = fidin.readlines()
            firstrow = otu_table_data[0].split('\t')
            sample_labels = firstrow[1:]
            sample_labels[len(sample_labels)-1] = sample_labels[len(sample_labels)-1].rstrip('\n')
            OTU_labels = [otu_table_data[i].split('\t')[0] for i in range(1,len(otu_table_data))]
            nOTUs = len(OTU_labels)
            nSamples = len(sample_labels)
            # Load OTU table row major order
            OTU_table_data = np.zeros((nOTUs, nSamples))
            for i in range(1,nOTUs+1):
                OTU_table_data[i-1,:] = otu_table_data[i].split('\t')[1:]

        # Initialize feature vectors and feature descriptions if they haven't already been initialized
        if not hasattr(self, 'feature_vectors'):
            self.feature_descriptions = {}
            self.feature_vectors = {}
            for sample_ID in sample_labels:
                self.feature_vectors[sample_ID] = []

        # Add OTU IDs to list of feature descriptions
        counter = len(self.feature_descriptions)
        for OTU_ID in OTU_labels:
            self.feature_descriptions[counter] = "OTU " + OTU_ID
            counter += 1
            
        # Input OTU abundances into feature vector
        for i in range(len(sample_labels)):
            sample_ID = sample_labels[i]
            try:
                self.feature_vectors[sample_ID] = self.feature_vectors[sample_ID] + OTU_table_data[:,i].tolist()
            except:
                print "No initialized feature vector found for sample ID " + sample_ID + ".  Do sample IDs differ between the OTU table and the metadata?"
                raise


    def LoadPredictedMetagenome(self, metapred_table_filename=0):
        # Load PiCRUST-predicted metagenome into feature vectors for each sample.
        if metapred_table_filename == 0:
            metapred_table_filename = self.metapred_table
        with open(metapred_table_filename,'r') as fidin:
            fileline = fidin.readlines()
            fileline = fileline[0]

            # Loop through the line to find KEGG pathway labels
            pathway_dictionary = {}
            start_str = '"rows": ['
            end_str = '],"columns":'
            starting_ind = fileline.find(start_str, 0)
            ending_ind = fileline.find(end_str,starting_ind)
            pathway_info = fileline[starting_ind + len(start_str) : ending_ind].split('},{')
            pathway_counter = 0
            for entry in pathway_info:
                entry_split = entry.split(',')
                pathway_name = entry_split[0].split()[1]
                pathway_dictionary[pathway_counter] = pathway_name
                pathway_counter += 1

            # Extract sample IDs
            sample_dictionary = {}
            start_str = '"columns": ['
            end_str = ']}'
            starting_ind = fileline.find(start_str, 0)
            ending_ind = fileline.find(end_str,starting_ind)
            sample_info = fileline[starting_ind + len(start_str) : ending_ind].split('},{')
            sample_counter = 0
            for entry in sample_info:
                entry_split = entry.split(',')
                sample_ID = entry_split[0].split(':')[1]
                sample_dictionary[sample_counter] = sample_ID.rstrip('"').lstrip(' "')
                sample_counter += 1

            # Extract pathway abundances for each sample
            feature_vectors = {}
            feature_names = {}
            starting_ind = fileline.find('[[', 0)
            starting_ind = starting_ind + 2
            ending_ind = fileline.find(']]', 0)
            metapred_data = fileline[starting_ind : ending_ind].split('],[')
            pathway_indices = []
            # init feature vectors
            for sample_index in sample_dictionary:
                sample_ID = sample_dictionary[sample_index]
                feature_vectors[sample_ID] = len(pathway_dictionary)*[0]
            # loop through data to extract pathway abundances
            for entry in metapred_data:
                entry_split = entry.split(',')
                pathway_index = entry_split[0]
                sample_index = entry_split[1]
                sample_ID = sample_dictionary[int(sample_index)]
                pathway_abundance = entry_split[2]
                feature_vectors[sample_ID][int(pathway_index)] = pathway_abundance

            # Load feature vectors into self
            for key in sample_dictionary:
                sample_ID = sample_dictionary[key]
                try:
                    self.feature_vectors[sample_ID] = self.feature_vectors[sample_ID] + feature_vectors[sample_ID]
                except:
                    print "No initialized feature vector found for sample ID " + sample_ID
                    raise

            # Add feature labels
            counter = len(self.feature_descriptions)
            for i in range(len(pathway_dictionary)):
                self.feature_descriptions[counter] = pathway_dictionary[i]
                counter += 1
            

    def LoadPhylogeneticFeatures(self, OTU_table_file=0, tree_file=0):
        # Load features from Phylogenetic tree routines
        if OTU_table_file == 0:
            OTU_table_file = self.otu_table_open
        if tree_file == 0:
            tree_file = self.phylogenetic_tree_open
        LminusR = phylo.ComputeTreeAbundances(OTU_table_file, tree_file)

        # Init feature arrays
        features = {}
        for sampleID in LminusR[LminusR.keys()[0]]:
            features[sampleID] = []

        # Add LminusR values as features and record the associated feature labels
        feature_counter = len(self.feature_descriptions)
        for node in LminusR:
            self.feature_descriptions[feature_counter] = "phylo_LminusR_" + node
            for sampleID in LminusR[node]:
                features[sampleID].append(LminusR[node][sampleID])
            feature_counter += 1

        # Append feature vectors to each sample feature vector
        for sampleID in features:
            self.feature_vectors[sampleID] = self.feature_vectors[sampleID] + features[sampleID]
            

    def LoadCollapsedTaxonomicFeatures(self, taxonomic_dict):
        # Adds the dictionary of collapsed taxonomies to the feature vector and descriptions
        # First, populate sampleIDs in feature vectors (if they don't exist yet)
        # Initialize feature vectors and feature descriptions if they haven't already been initialized
        sampleIDs = taxonomic_dict[taxonomic_dict.keys()[0]].keys()
        if not hasattr(self, 'feature_vectors'):
            self.feature_descriptions = {}
            self.feature_vectors = {}
            for sample_ID in sampleIDs:
                self.feature_vectors[sample_ID] = []
        feature_counter = len(self.feature_descriptions)

        for taxon in taxonomic_dict:
            self.feature_descriptions[feature_counter] = taxon
            for sampleID in taxonomic_dict[taxon]:
                self.feature_vectors[sampleID].append(taxonomic_dict[taxon][sampleID])
            feature_counter += 1        

