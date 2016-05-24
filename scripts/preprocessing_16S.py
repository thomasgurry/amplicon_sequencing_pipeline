"""

OVERVIEW: 

Python module for pre-processing FASTQ data into OTU tables.

"""


import numpy as np
import sys
import os, sys
import util
import Formatting
from bidict import *
import pandas as pd
from collections import defaultdict

def length_stats_fastq(fastq_in):
    # Returns full sequence length, and  5th percentile of read length for a 100000 sample from a FASTQ file.
    iter_seq = util.iter_fsq
    x = []
    counter = 0
    for record in iter_seq(fastq_in):
        sid, seq = record[:2]
        counter = counter + 1
        if(counter > 100000):
            break
        x.append(len(seq))
    x = np.array(x)
    return [np.amax(x), np.percentile(x,5)]

def trim_quality((fastq_in, fastq_out, ascii_encoding, quality_trim)):
    # Trims to desired quality level
    str1 = '/home/ubuntu/bin/usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(quality_trim) + ' -fastq_ascii ' + str(ascii_encoding) + ' -fastqout ' + fastq_out
    os.system(str1)
    return None


def trim_quality_deprecated((fastq_in, fastq_out, ascii_encoding)):
    # Finds maximum Q-score cut-off where 95% of reads are over 200 bases, if possible.  If not, selects Q=5 and returns a warning.
    Qvals = range(5,11)
    bestQ = 0

    for Q in Qvals:
        # Trim to Q
        print "[[ Quality trimming ]] Quality trimming with Q=" + str(Q)        
        str1 = '/home/ubuntu/bin/usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(Q) + ' -fastq_ascii ' + str(ascii_encoding) + ' -fastqout ' + fastq_out
        os.system(str1)
        # Check length distribution
        try:
            [full_length, fifthPercentile] = length_stats_fastq(fastq_out)
            print "[[ Quality trimming ]] 5th percentile length: " + str(fifthPercentile)
            if(fifthPercentile >= 0.8*float(full_length)):
                bestQ = Q
                bestFifthPercentile = fifthPercentile
            else:
                continue
        except:
            continue
    if (bestQ == 0):
        print "[[ Quality trimming ]] ERROR!!  Could not obtain 95% of reads over 200 base pairs with quality score cut-off of at least 5.  Check sequencing quality of the data.  Proceeding with Q=5..."
        str1 = '/home/ubuntu/bin/usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(5) + ' -fastq_ascii ' + str(ascii_encoding) + ' -fastqout ' + fastq_out
        os.system(str1)
        [new_full_length, fifthPercentile] = length_stats_fastq(fastq_out)
        print "[[ Quality trimming ]] Input file: " + fastq_in
        print "[[ Quality trimming ]] ASCII encoding used: " + str(ascii_encoding)
        print "[[ Quality trimming ]] Trimmed sequences with quality cut-off of Q = " + str(5) + "."  
        print "[[ Quality trimming ]] 5th percentile of read length: " + str(fifthPercentile) + "."
        print "[[ Quality trimming ]] Trimmed file: " + fastq_out
        print "[[ Quality trimming ]] Complete."
        return None
           
    else:
        str1 = '/home/ubuntu/bin/usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(bestQ) + ' -fastq_ascii ' + str(ascii_encoding) + ' -fastqout ' + fastq_out
        os.system(str1)
        print "[[ Quality trimming ]] Input file: " + fastq_in
        print "[[ Quality trimming ]] ASCII encoding used: " + str(ascii_encoding)
        print "[[ Quality trimming ]] Trimmed sequences with quality cut-off of Q = " + str(bestQ) + "."  
        print "[[ Quality trimming ]] 5th percentile of read length: " + str(bestFifthPercentile) + "."
        print "[[ Quality trimming ]] Trimmed file: " + fastq_out
        print "[[ Quality trimming ]] Complete."
        return None

def trim_length_fasta((fasta_in, fasta_out, length)):
    # Trims FASTA files to uniform length
    str1 = '/home/ubuntu/bin/usearch8 -fastx_truncate ' + fasta_in + ' -trunclen ' + str(length) + ' -fastaout ' + fasta_out
    os.system(str1)

def trim_length_fastq((fastq_in, fastq_out, length, ascii_encoding)):
    # Trims FASTQ files to usearch -fastq_chars obio.raw.fsq uniform length and filters by maximum expected error
    # Takes as input the ascii encoding (33 or 64 currently supported)
    str1 = '/home/ubuntu/bin/usearch8 -fastq_filter ' + fastq_in + ' -fastq_trunclen ' + str(length) + ' -fastq_ascii ' + str(ascii_encoding) + ' -fastqout ' + fastq_out
    os.system(str1)
    # Check for how many reads were thrown out
    statinfoIN = os.stat(fastq_in)
    input_filesize = float(statinfoIN.st_size)
    statinfoOUT = os.stat(fastq_out)
    output_filesize = float(statinfoOUT.st_size)
    percent_thrown_out = 100*(output_filesize / input_filesize)
    print "[[ Length trimming ]] Input file: " + fastq_in
    print "[[ Length trimming ]] Trimmed sequences with length less than " + str(length) + " and maximum expected error of 0.25."
    print "[[ Length trimming ]] Threw out " + str(percent_thrown_out) + " % of reads."
    print "[[ Length trimming ]] Complete."
    return None
 
def remove_primers((fastq_in, fastq_out, primers_file)):
    # Remove primers from FASTQ file
    print "[[ Primer trimming ]] ..."
    os.system('python ~/scripts/1.remove_primers.py -q ' + fastq_in + ' -l ' + primers_file + ' -d 1 -o ' + fastq_out)
    print "[[ Primer trimming ]] Complete."
    return None


def split_by_barcodes((fastq_in, fastq_out, barcodes_map, mode)):
    # Split by barcodes
    print "[[ Splitting by barcodes ]] ..."
    os.system('python ~/scripts/2.split_by_barcodes.py -q ' + fastq_in + ' -b ' + barcodes_map + ' -B tab -d 1 --mode ' + mode + ' -o ' + fastq_out)
    return None

def split_by_barcodes_FASTQ((fastq_in, fastq_out, barcodes_map, mode)):
    # Split by barcodes
    print "[[ Splitting by barcodes ]] ..."
    os.system('python ~/scripts/2.split_by_barcodes.py -q ' + fastq_in + ' -b ' + barcodes_map + ' -B tab -d 1 --mode ' + mode + ' -o ' + fastq_out)
    return None

def split_by_barcodes_FASTA((fasta_in, fasta_out, barcodes_map, mode)):
    # Split by barcodes
    print "[[ Splitting by barcodes ]] ..."
    os.system('python ~/scripts/2.split_by_barcodes.py -f ' + fasta_in + ' -b ' + barcodes_map + ' -B tab -d 1 --mode ' + mode + ' -o ' + fasta_out)
    return None


def replace_seqIDs_for_demultiplexed_files((fastq_in, fastq_out, sampleID)):
    # Relabels all seqIDs with the provided sample ID.  For use when a single raw file has reads only for one sample and this is known.
    outfile = open(fastq_out, 'w')
    iter_fsq = util.iter_fsq
    seq_counter = 0
    for record in iter_fsq(fastq_in):
        sid = record[0][1:] # id
        seq = record[1] # sequence
        record[0] = '@' + sampleID + '_' + str(seq_counter)
        outfile.write('\n'.join(record) + '\n')
        seq_counter += 1
    outfile.close()


def dereplicate_and_sort(fasta_in, fasta_out, OTU_database, separator, processing_summary_file, min_count):
    # Dereplicate and sort sequences by size
    print "[[ Dereplicating and sorting ]] Discarding sequences with fewer than " + str(min_count) + " reads"
    os.system('python ~/scripts/3.dereplicate.py -f ' + fasta_in + " -s '" + separator + "' -o " + OTU_database + ' -d ' + fasta_out + ' -P ' + processing_summary_file + ' -M ' + str(min_count))
    print "[[ Dereplicating and sorting ]] Complete."
    return None


def remove_chimeras_and_cluster_OTUs(fasta_in, OTU_sequences_fasta, clustering_results, cluster_percentage=97.0, relabel=False):
    # Remove chimeric sequences and then cluster OTUs with default similarity of 97%.
    print "[[ Removing chimeras and clustering OTUs ]] ..."
    max_cluster_diff = 100.0 - float(cluster_percentage)
    if relabel == True:
        os.system('/home/ubuntu/bin/usearch8 -cluster_otus ' + fasta_in + ' -otus ' + OTU_sequences_fasta + ' -otu_radius_pct ' + str(max_cluster_diff) + ' -sizein -uparseout ' + clustering_results + ' -relabel denovo')
    else:
        os.system('/home/ubuntu/bin/usearch8 -cluster_otus ' + fasta_in + ' -otus ' + OTU_sequences_fasta + ' -otu_radius_pct ' + str(max_cluster_diff) + ' -sizein -uparseout ' + clustering_results)

    print "[[ Removing chimeras and clustering OTUs ]] Complete."
    return None


def separate_GG_reads(fasta_dereplicated, OTU_GG_dict, output_GG_reads, output_denovo_reads):
    # Takes as input a dereplicated set of reads and an OTU alignment dictionary, and separates the reads into GG-referenced and not GG-referenced
    iter_fst = util.iter_fst
    with open(output_GG_reads, 'w') as GGfid:
        with open(output_denovo_reads, 'w') as dnfid:
            for record in iter_fst(fasta_dereplicated):
                sid = record[0][1:] # sequence ID
                if sid in OTU_GG_dict:
                    GGfid.write('\n'.join(record) + '\n')
                else:
                    newsid = sid
                    record[0] = record[0][0] + newsid
                    dnfid.write('\n'.join(record) + '\n')


def collapse_oligotypes(oligotype_table_filename, output_OTU_table):
    # Collapse oligotypes onto OTUs (1.1, 1.2, 1.3 --> 1) and print an OTU table
    x = pd.read_csv(oligotype_table_filename, sep='\t')
    oligos = x['#Oligotype']
    columns = x.columns

    oligo_dict = {}
    OTU_dict = {}

    for i in range(len(oligos)):
        counts = np.array(x.ix[i])
        oligo_dict[oligos[i]] = counts[1:]
        OTU_ID = str(str(oligos[i]).split('.')[0])
        if OTU_ID not in OTU_dict:
            OTU_dict[OTU_ID] = counts[1:]
        else:
            OTU_dict[OTU_ID] = OTU_dict[OTU_ID] + counts[1:]

    # Write OTU table
    with open(output_OTU_table, 'w') as fid:
        firstline = 'OTU_ID' + '\t' + '\t'.join(columns[1:])
        fid.write(firstline + '\n')
        for OTU in OTU_dict:
            abundance_string = '\t'.join(OTU_dict[OTU].astype(int).astype(str))
            line = OTU + '\t' + abundance_string
            fid.write(line + '\n')
    

def build_GG_OTU_table(dereplication_map, OTU_GG_dict, OTU_table_gg):
    # Builds a GG-referenced OTU table from a list of sequences (FASTA), a dereplication map (from dereplicate_and_sort()) and a dictionary mapping the fasta sequence IDs to OTU GG IDs.
    gg_dereplication_map = dereplication_map + '.gg'
    OTU_GG_dict_new = {}
    for key in OTU_GG_dict:
        OTU_GG_dict_new[key.split(';')[0]] = OTU_GG_dict[key]
    
    # write a new dereplication map which has only the seqs that mapped to GGID
    # also store this info in gg_derep_dict, which has {sampleID: {original_OTU: count}}
    gg_derep_dict = defaultdict(lambda: defaultdict(str))
    with open(gg_dereplication_map,'w') as newfid:
        with open(dereplication_map,'r') as fid:
            derepmap = fid.readlines()
            new_derepmap = [line for line in derepmap if line.split()[0] in OTU_GG_dict_new]
            for line in new_derepmap:
                linespl = line.split('\t')
                # First, grab the original ID and sample map info to update the gg_derep_dict (bc dicts can't handle having duplicate keys)
                seqid = linespl[0]
                smpl_map = linespl[1].split(' ')
                smpl_map = {s.split(':')[0]: s.split(':')[1] for s in smpl_map}
                for sid in smpl_map:
                    gg_derep_dict[sid][seqid] = smpl_map[sid]
                # Now update the line with new OTU GG ID and write to gg derep file
                linespl[0] = OTU_GG_dict_new[linespl[0]]
                linejoined = '\t'.join(linespl)
                newfid.write(linejoined)

    # write GG OTU table
    with open(OTU_table_gg, 'w') as f:
        gg_otus = list(set(OTU_GG_dict.values()))
        f.write('sample\t' + '\t'.join(gg_otus) + '\n')
        for smpl in gg_derep_dict:
            f.write(smpl + '\t')
            seq_dict = dict.fromkeys(gg_otus, 0)
            for orig_otu in gg_derep_dict[smpl]:
                new_otu = OTU_GG_dict_new[orig_otu]
                seq_dict[new_otu] = seq_dict[new_otu] + int(gg_derep_dict[smpl][orig_otu])
            counts = [str(seq_dict[otu]) for otu in gg_otus]
            f.write('\t'.join(counts) + '\n')


def parse_multihit_alignment(uc_file):
    # Parses a UC file from a usearch alignment where -maxaccepts is set to 10, i.e. there are a maximum of ten hits per query.
    with open(uc_file, 'r') as fid:
        all_lines = fid.readlines()
    alignment_dict = {}
    for line in all_lines:
        query = line.split()[8]
        if query not in alignment_dict.keys():
            alignment_dict[query] = []
        hit = line.split()[9]
        alignment_dict[query].append(hit)
    return alignment_dict


def parse_alignment(alignment_file):
    # Parses an alignment of OTU sequences against a reference database.  Returns a dict of query indices (OTU IDs) and associated reference database IDs: {'OTU_ID': DB_ID}
    with open(alignment_file, 'r') as fid:
        all_lines = fid.readlines()
        line_nums = [i for i in range(len(all_lines)) if all_lines[i][:6] == " Query"]
        alignment_dict = {}
        for line_num in line_nums:
            line1 = all_lines[line_num]
            line2 = all_lines[line_num + 1]
            query = line1.split()[2]
            query = query[1:]
            target = line2.split()[2]
            target = target[1:]
            alignment_dict[query] = target
    return alignment_dict

def renumber_sequences(fasta_files, separator):
    # Renumbers sequences IDs for each sample so sampleID_1, sampleID_2 etc. only occur once
    def update_numbers(fasta_in, fasta_out, sample_counts):
        iter_fst = util.iter_fst
        with open(fasta_out, 'w') as outfile:
            for record in iter_fst(fasta_in):
                sid = record[0][1:] # id
                sid = sid.split(separator)
                sampleID = ''.join(sid[:len(sid)-1])
                seq = record[1] # sequence
                sample_count = sid[len(sid)-1]
                sample_counts[sampleID] = sample_counts.get(sampleID, 0) + 1 # Increment sample count
                new_sid = '>%s_%d' %(sampleID, sample_counts[sampleID])
                record[0] = new_sid
                outfile.write('\n'.join(record) + '\n')
        return sample_counts

    print "[[ Renumbering sequences ]] ..."
    sample_counts = {}
    for filename in fasta_files:
        sample_counts = update_numbers(filename, filename + '.tmp', sample_counts)
        os.system('mv ' + filename + '.tmp ' + filename)


def pull_counts(raw_sequence_filename, separator):
    # Builds a dict of counts of each unique sequence in each sample, x[seq][sampleID]
    x = {}
    fn = raw_sequence_filename
    iter_fst = util.iter_fst
    samples = []
    for record in iter_fst(fn):
        [sid, seq] = record[:2]
        sid = sid[1:]
        sa = sid.split(separator)
        sa = separator.join(sa[:len(sa)-1])
        if seq not in x:
            x[seq] = {}
        if sa not in x[seq]:
            x[seq][sa] = 0
        if sa not in samples:
            samples.append(sa)
        x[seq][sa] += 1
    return x, samples 


def compute_oligotype_table(raw_trimmed, raw_dereplicated, clustering_file, separator, oligotype_table_filename):
    # Inputs: 
    #       'raw_trimmed' = raw reads, trimmed to final form (FASTA)
    #       'raw_dereplicated' = dereplicated reads (FASTA)
    #       'clustering_file' = output clustering file from 'usearch8 -cluster_otus' (-parseout option) 
    #       'separator' = separator character, e.g. '_' for '>sampleID_sequenceNumber' sequence IDs
    #       'oligotype_table_filename' = output filename
    #
    # 1. get sequence counts
    x, samples = pull_counts(raw_trimmed, separator)

    # 2. populate sequence lookup 
    sequence_lookup = bidict({})    # seq <--> seqID
    iter_fst = util.iter_fst
    for record in iter_fst(raw_dereplicated):
        [sid, seq] = record[:2]
        sid = int(sid[1:].split(';')[0])
        sequence_lookup[seq] = sid

    # 3.  Populate clustering_lookup (from otu_clustering.tab)
    clustering_lookup = {}  # seqID <--> 'otu' or 'match'
    OTU_lookup = {}         # seqID <--> OTU_ID centroid
    with open(clustering_file, 'r') as fid:
        all_lines = fid.readlines()
        for line in all_lines:
            split_line = line.split()
            seqID = int(split_line[0].split(';')[0])
            clustering_lookup[seqID] = split_line[1]
            if split_line[1] == 'match' or split_line[1] == "otu":
                OTU_lookup[seqID] = split_line[4]
    
    # 4.  Populate dictionaries with each sequence within an OTU.  Each of the three dictionaries contain lists whose entries are ordered in the same manner in each dict.

    OTU_oligos = {}         # OTU_ID <--> ['ACAGT','ACAAT', ...] 
    OTU_original_seqIDs = {}  # OTU_ID <--> [seqID1, seqID2, ...] (original sequence IDs for oligotypes)
    OTU_oligo_IDs = {}      # OTU_ID <--> [0, 1, ...] (oligotype IDs)
    for seq in sequence_lookup:
        seqID = sequence_lookup[seq]
        if clustering_lookup[seqID] != "chimera":
            OTU_centroid = OTU_lookup[seqID]
            if OTU_centroid not in OTU_oligos:
                OTU_oligos[OTU_centroid] = []
                OTU_oligo_IDs[OTU_centroid] = []
                OTU_original_seqIDs[OTU_centroid] = []
            if seq not in OTU_oligos[OTU_centroid]:
                OTU_oligos[OTU_centroid].append(seq)
                OTU_original_seqIDs[OTU_centroid].append(seqID)
                if len(OTU_oligo_IDs[OTU_centroid]) > 0:
                    OTU_oligo_IDs[OTU_centroid].append(OTU_oligo_IDs[OTU_centroid][-1] + 1)
                else:
                    OTU_oligo_IDs[OTU_centroid].append(0)


    # Create a full list of oligotypes
    oligotype_list = []
    for OTU in OTU_oligo_IDs:
        for oligoID in OTU_oligo_IDs[OTU]:
            full_oligotype_ID = str(OTU) + '.' + str(oligoID)
            oligotype_list.append(full_oligotype_ID)

    # Get counts for each oligotype
    oligotype_counts = {}
    # Loop through each sequence in x and for each OTU,
    for seq in x:
        # assign oligotype_ID
        try:
            seqID = sequence_lookup[seq]
            if clustering_lookup[seqID] != "chimera":
                OTU_centroid = OTU_lookup[seqID]
                # Look up which oligotype ID it is
                index = OTU_oligos[OTU_centroid].index(seq)
                oligo_ID = OTU_oligo_IDs[OTU_centroid][index]
                full_oligotype_ID = str(OTU_centroid) + '.' + str(oligo_ID)
                # Get counts for each sample
                oligotype_counts[full_oligotype_ID] = x[seq]
        except:
            continue


    # Create an oligotype table, where each row is an oligotype (e.g. 1.0 or 1.1) and each column is a sample    
    with open(oligotype_table_filename, 'w') as fid:
        firstline = "#Oligotype" + '\t' + '\t'.join(samples)
        fid.write(firstline + '\n')
        for full_oligotype_ID in oligotype_counts:
            oligo_abundances_per_sample = ['0']*len(samples)
            for sampleID in oligotype_counts[full_oligotype_ID]:
                index = samples.index(sampleID) # get index in samples vector
                oligo_abundances_per_sample[index] = str(oligotype_counts[full_oligotype_ID][sampleID])
            line = full_oligotype_ID + '\t' + '\t'.join(oligo_abundances_per_sample)
            fid.write(line + '\n')


def concatenate_OTU_tables(OTU_table_1, OTU_table_2, combined_OTU_table):
    # Concatenates two classic format OTU tables (with the same samples as columns)
    x = pd.read_csv(OTU_table_1, sep='\t')
    y = pd.read_csv(OTU_table_2, sep='\t')
    z = pd.concat([x, y])
    cols = z.columns.tolist()
    newcols = ['OTU_ID']
    for col in cols:
        if col != 'OTU_ID':
            newcols.append(col)
    z.to_csv(combined_OTU_table, columns=newcols, sep='\t', index=False)


def build_OTU_table_from_alignments(dereplication_map, OTU_dict, OTU_table):
    # Builds a database-referenced OTU table for a list of sequences (FASTA), a dereplication map (from dereplicate_and_sort()) and a dictionary mapping the fasta sequence IDs to the database IDs.
    new_dereplication_map = dereplication_map + '.new'
    OTU_dict_new = {}
    for key in OTU_dict:
        OTU_dict_new[key.split(';')[0]] = OTU_dict[key]
    
    # write a new dereplication map which has only the seqs that mapped to ITS
    # also store this info in derep_dict, which has {sampleID: {original_OTU: count}}
    derep_dict = defaultdict(lambda: defaultdict(str))
    with open(new_dereplication_map,'w') as newfid:
        with open(dereplication_map,'r') as fid:
            derepmap = fid.readlines()
            new_derepmap = [line for line in derepmap if line.split()[0] in OTU_dict_new]
            for line in new_derepmap:
                linespl = line.split('\t')
                # First, grab the original ID and sample map info to update the derep_dict (bc dicts can't handle having duplicate keys)
                seqid = linespl[0]
                smpl_map = linespl[1].split(' ')
                smpl_map = {s.split(':')[0]: s.split(':')[1] for s in smpl_map}
                for sid in smpl_map:
                    derep_dict[sid][seqid] = smpl_map[sid]
                # Now update the line with new OTU ID and write to derep file
                linespl[0] = OTU_dict_new[linespl[0]]
                linejoined = '\t'.join(linespl)
                newfid.write(linejoined)

    # write OTU table
    with open(OTU_table, 'w') as f:
        DB_otus = list(set(OTU_dict.values()))
        f.write('sample\t' + '\t'.join(DB_otus) + '\n')
        for smpl in derep_dict:
            f.write(smpl + '\t')
            seq_dict = dict.fromkeys(DB_otus, 0)
            for orig_otu in derep_dict[smpl]:
                new_otu = OTU_dict_new[orig_otu]
                seq_dict[new_otu] = seq_dict[new_otu] + int(derep_dict[smpl][orig_otu])
            counts = [str(seq_dict[otu]) for otu in DB_otus]
            f.write('\t'.join(counts) + '\n')


def collapse_alignment_dict(alignment_dict, minimum_hits, db='GG'):
    # Assign each sequence ID to a single consensus alignment based on a minimum number of top hits being in agreement.
    # Assumes format of GreenGenes database by default.  Options: db=['GG','UNITE']
    taxonomic_levels  = ['k','p','c','o','f','g','s','d']
    new_dict = {}
    for key in alignment_dict:
        assignments = alignment_dict[key]
        seqid = key.split(';')[0]
        if len(assignments) >= minimum_hits and assignments[0] != '*':
            if db == 'UNITE':
                assignments = [ass.split('|')[1] for ass in assignments]            
            consensus_assignment = []
            # Check each level to build consensus
            unidentified = 0
            for i in range(7):
                try:
                    tax_level_assignments = [ass.split(';')[i] for ass in assignments]
                    if len(np.unique(tax_level_assignments[:minimum_hits])) == 1 and unidentified == 0:
                        consensus_assignment.append(tax_level_assignments[0])
                    else:
                        consensus_assignment.append(taxonomic_levels[i] + '__')
                        unidentified = 1

                except:
                    consensus_assignment.append(taxonomic_levels[i] + '__')
            consensus_assignment.append('d__derep' + seqid)
            new_dict[key] = ';'.join(consensus_assignment)
        else:
            new_dict[key] = 'k__;p__;c__;o__;f__;g__;s__;d__derep' + seqid
    return new_dict


def RDP_classify(OTU_sequences_fasta, RDP_classifications_file, amplicon_type='16S'):
    # Function to run Scott's script for RDP classification.
    if amplicon_type == '16S':
        cmd_str = 'python ~/scripts/rdp_classify.py ' + OTU_sequences_fasta + ' ' + RDP_classifications_file + ' --gene 16srrna'
    elif amplicon_type == 'ITS':
        cmd_str = 'python ~/scripts/rdp_classify.py ' + OTU_sequences_fasta + ' ' + RDP_classifications_file + ' --gene fungalits_unite'
    os.system(cmd_str)
    return None


def parse_RDP_classifications(rdp_classifications, cutoff):
    # Function to parse RDP classifications based on a desired probability cutoff
    with open(rdp_classifications,'r') as fid:
        all_lines = fid.readlines()
    RDP_assignments = {}

    for line in all_lines:
        line = line.strip().split('\t')
        seqID = line[0]
        assignment = ''
        # Domain
        try:
            i = line.index('domain')
            domain = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += 'k__' + domain
            else:
                assignment += 'k__'
        except:
            assignment += 'k__'
        # Phylum
        try:
            i = line.index('phylum')
            phylum = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';p__' + phylum
            else:
                assignment += ';p__'
        except:
            assignment += ';p__'
        # Class
        try:
            i = line.index('class')
            class_assignment = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';c__' + class_assignment
            else:
                assignment += ';c__'
        except:
            assignment += ';c__'
        # Order
        try:
            i = line.index('order')
            order = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';o__' + order
            else:
                assignment += ';o__'
        except:
            assignment += ';o__'
        # Family
        try:
            i = line.index('family')
            family = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';f__' + family
            else:
                assignment += ';f__'
        except:
            assignment += ';f__'
        # Genus
        try:
            i = line.index('genus')
            genus = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';g__' + genus
            else:
                assignment += ';g__'
        except:
            assignment += ';g__'
        # Species
        try:
            i = line.index('species')
            species = line[i - 1]
            prob = float(line[i + 1])
            if prob >= cutoff:
                assignment += ';s__' + species
            else:
                assignment += ';s__'
        except:
            assignment += ';s__'
        # Add the seqID as a unique identifier (i.e. keep the 'denovo' OTU ID)
        assignment += ';d__' + str(seqID)
        # Remove all quotes in the assignment
        RDP_assignments[seqID] = assignment.replace('"','')

    return RDP_assignments


def relabel_denovo_OTUs_with_RDP(OTU_table_denovo, RDP_assignments):
    # Relabels the OTU IDs with their latin names
    otu_table = pd.read_csv(OTU_table_denovo,sep='\t')
    OTU_IDs = otu_table['OTU_ID'].tolist()
    for i in range(len(OTU_IDs)):
        OTU_IDs[i] = RDP_assignments[OTU_IDs[i]]
    otu_table['OTU_ID'] = OTU_IDs
    otu_table.to_csv(OTU_table_denovo + '.rdp_assigned', sep='\t', index=False)
