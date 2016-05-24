"""

OVERVIEW: 

Python module for depositing processed data into The Database.
Takes a summary file as input.

"""

import os, sys
import numpy as np
from optparse import OptionParser
import csv
from SummaryParser import *
import MySQLdb
import util

class DepositDB():
    def __init__(self, data_directory):
        self.summary_file = data_directory + '/summary_file.txt'
        #self.DepositOTUTable()

    def DepositOTUSequences(self):
        data_summary = SummaryParser(self.summary_file)
        data_summary.ReadSummaryFile()        
        # Load sequences
        try:
            otu_sequence_fasta_filename = data_summary.attribute_value_16S['OTU_SEQUENCES_FASTA']
            iter_fst = util.iter_fst
            otu_seqs = []
            otu_IDs = []
            for record in iter_fst(otu_sequence_fasta_filename):
                [sid, seq] = record[:2]
                sid = sid[1:]
                otu_seqs.append(seq)
                otu_IDs.append(sid)        
        except:
            otu_sequence_table_filename = data_summary.attribute_value_16S['OTU_SEQUENCES_TABLE'] 
            otu_sequence_table_data = open(otu_sequence_table_filename, 'r').readlines()
            otu_IDs = [otuinfo.split('\t')[0] for otuinfo in otu_sequence_table_data]
            otu_IDs = otu_IDs[1:]
            otu_seqs = [otuinfo.split('\t')[1].rstrip('\n') for otuinfo in otu_sequence_table_data]
            otu_seqs = otu_seqs[1:]

        ##########################
        # OTU_sequences database #
        ##########################
        db = MySQLdb.connect(host="mbitdb1.cwnerb5rh4qg.us-east-1.rds.amazonaws.com", # your host, usually localhost
                     user="cmitaws", # your username
                     passwd="microbiome", # your password
                     db="OTU_sequences") # name of the data base
        insert_OTU_seq = ("INSERT INTO turnbaugh_mouse_diet_2015 "
                          "(OTU_ID, OTU_sequence, latin_name)"
                          "VALUES (%(OTU_ID)s, %(OTU_sequence)s, %(latin_name)s)"
        )
        cursor = db.cursor()

        # Create OTU table in database 'OTU_tables'
        mysql_cmd = "create table " + data_summary.datasetID + " (OTU_ID VARCHAR(20), OTU_sequence VARCHAR(500), latin_name VARCHAR(300));"
        cursor.execute(mysql_cmd)

        # Deposit OTU sequences
        for i in range(len(otu_seqs)):
            # Do stuff here
            data_OTU_seq = {
                'OTU_ID': otu_IDs[i],
                'OTU_sequence': otu_seqs[i],
                'latin_name': "NULL"
            }
            cursor.execute(insert_OTU_seq, data_OTU_seq)
        db.commit()


    def DepositOTUTable(self):
        data_summary = SummaryParser(self.summary_file)
        data_summary.ReadSummaryFile()

        # Extract OTU counts
        otu_table_filename = data_summary.attribute_value_16S['OTU_TABLE']
        with open(otu_table_filename, 'r') as fid:
            otu_table_data = fid.readlines()
            first_sample = otu_table_data[1].split('\t')
            nOTUs = len(first_sample) - 1

            sample_labels = [otu_table_data[i].split('\t')[0] for i in range(1,len(otu_table_data))]
            nSamples = len(sample_labels)-1

        # Load OTU table
        OTU_table = np.loadtxt(otu_table_filename, skiprows=1, usecols=range(1,nOTUs+1))
        
        with open(otu_table_filename, 'r') as fid:
            otu_table_data = fid.readlines()
            OTU_labels = otu_table_data[0].split('\t')
            OTU_labels = OTU_labels[1:]
            OTU_labels[len(OTU_labels)-1] = OTU_labels[len(OTU_labels)-1].rstrip('\n')
        
        # Normalize OTU counts
        row_sums = OTU_table.sum(axis=1)
        OTU_table_norm = OTU_table / row_sums[:, np.newaxis]

        # Check if table has already been deposited.  If not, create an entry.
        #
        # -need sample_info = studyID [sampleIDs / disease_labels / keywords] database
        # -need study_info = studyIDs [studyID / study_name / study_description / keywords] database
        # -need OTU_sequences = studyID [OTU_ID / latin_name / sequence] database
        # 


        ###################### 
        # OTU_table database #
        ######################
        db = MySQLdb.connect(host="mbitdb1.cwnerb5rh4qg.us-east-1.rds.amazonaws.com", # your host, usually localhost
                     user="cmitaws", # your username
                     passwd="microbiome", # your password
                     db="OTU_tables") # name of the data base
        insert_OTU = ("INSERT INTO CRC_Zhao_2012 " 
                      "(sample_ID, OTU_ID, abundance_counts, abundance_norm) "
                      "VALUES (%(sample_ID)s, %(OTU_ID)s, %(abundance_counts)s, %(abundance_norm)s)")
        cursor = db.cursor()

        # Create OTU table in database 'OTU_tables'
        mysql_cmd = "create table " + data_summary.datasetID + " (sample_ID VARCHAR(20), OTU_ID VARCHAR(20), abundance_counts INT(11), abundance_norm float);"
        cursor.execute(mysql_cmd)
        
        for i in range(nSamples):
            for j in range(nOTUs):
                # Do stuff here
                data_OTU = {
                    'studyID': data_summary.datasetID,
                    'sample_ID': sample_labels[i],
                    'OTU_ID': OTU_labels[j], 
                    'abundance_counts': OTU_table[i,j], 
                    'abundance_norm': OTU_table_norm[i,j],
                }
                cursor.execute(insert_OTU, data_OTU)
        db.commit()


    def DepositSampleInfo(self):

        ######################## 
        # Sample_info database #
        ########################

        # Insert into DB
        db = MySQLdb.connect(host="mbitdb1.cwnerb5rh4qg.us-east-1.rds.amazonaws.com", # your host, usually localhost
                     user="cmitaws", # your username
                     passwd="microbiome", # your password
                     db="Sample_info") # name of the data base
        insert_OTU = ("INSERT INTO turnbaugh_mouse_diet_2015 "
                      "(sample_ID, disease_labels, keywords)" 
                      "VALUES (%(sample_ID)s, %(disease_labels)s, %(keywords)s)")
        cursor = db.cursor()


        # CREATE sample info table
        mysql_cmd1 = "create table " + data_summary.datasetID + " (sampleID VARCHAR(20), disease_labels VARCHAR(30), keywords VARCHAR(200));"

        # Deposit sample info
        
