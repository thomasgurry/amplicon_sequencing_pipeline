"""

OVERVIEW: 

Python module for extracting taxonomic information from closed-reference OTU tables in classic dense format.

"""

import os, sys
import numpy as np
import Formatting
import pandas as pd

def GGID_lookup(GG_ID, taxonomy_file=0):
    # Takes as input a GreenGenes ID and returns a latin name.  Optionally, can specify which GG taxonomy file to use.  By default uses 97% reference set.
    if(taxonomy_file == 0):
        taxonomy_file = '/home/ubuntu/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt'

    # Extract latin names for OTU IDs
    with open(taxonomy_file, 'r') as fid:
        all_lines = fid.readlines()
        gg_OTU_IDs = [line.split()[0] for line in all_lines]
    
    ID_index = gg_OTU_IDs.index(GG_ID)
    return "".join(all_lines[ID_index].split()[1:])
        

def convert_GGIDs_to_latin_names(OTU_table, new_OTU_table, taxonomy_file=0):
    # Takes as input a closed-reference OTU table with GreenGenes IDs and returns the same table but with latin names instead of GG IDs.
    if(taxonomy_file == 0):
        taxonomy_file = '/home/ubuntu/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt'

    # Extract OTU IDs
    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]

    # Extract latin names for OTU IDs
    with open(taxonomy_file, 'r') as fid:
        all_lines = fid.readlines()
        gg_OTU_IDs = [line.split()[0] for line in all_lines]
        indices_of_OTUs_in_gg = [gg_OTU_IDs.index(OTU) for OTU in OTU_IDs]
        OTU_latin_names = ["".join(all_lines[i].split()[1:]) for i in indices_of_OTUs_in_gg]

    # Rename OTUs by latin names in OTU table
    for i in range(1,len(OTU_table_lines)):
        line = OTU_table_lines[i].split()
        line[0] = OTU_latin_names[i-1]
        OTU_table_lines[i] = '\t'.join(line)

    # Rewrite OTU table
    with open(new_OTU_table, 'w') as fid:
        fid.write(OTU_table_lines[0])
        for i in range(1,len(OTU_table_lines)):
            fid.write(OTU_table_lines[i] + '\n')

def collapse_taxonomic_contents(OTU_table, taxonomic_level, parent_node=0):
    # Takes as input an OTU table with latin names as OTU IDs,
    # the taxonomic level to which all taxonomies should be collapsed.
    # Accepts: kingdom, phylum, class, order, family, genus, species
    # e.g.  collapse_taxonomic_contents(OTU_table, class, Firmicutes)
    # If no parent node is specified, returns all clades of the given taxonomic level.

    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]

    # Collapse to the right level        
    if(taxonomic_level == "domain"):
        OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "phylum"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[0][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[0]) > 3]
    if(taxonomic_level == "class"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[1][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[1]) > 3]
    if(taxonomic_level == "order"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[2][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[2]) > 3]
    if(taxonomic_level == "family"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[3][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[3]) > 3]
    if(taxonomic_level == "genus"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[4][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[4]) > 3]
    if(taxonomic_level == "species"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[5][3:] for OTU_ID in OTU_IDs if len(OTU_ID.split(';')[5]) > 3]
        

    # Get indices of each unique taxon
    taxa_indices = {}
    if parent_node != 0:
        for i in range(len(OTU_taxa)):
            if (OTU_taxa[i] not in taxa_indices) and (OTU_taxa_parents[i] == parent_node) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
                taxa_indices[OTU_taxa[i]] = []
            if (OTU_taxa_parents[i] == parent_node) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
                taxa_indices[OTU_taxa[i]].append(i)            
    elif parent_node == 0:
        for i in range(len(OTU_taxa)):
            if (OTU_taxa[i] not in taxa_indices) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
                taxa_indices[OTU_taxa[i]] = []
            if (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
                taxa_indices[OTU_taxa[i]].append(i)

    # Read in OTU abundances
    sample_labels = OTU_table_lines[0].split()[1:]
    nOTUs = len(OTU_IDs)
    nSamples = len(sample_labels)

    # Load OTU table row major order
    OTU_table_counts = np.zeros((nOTUs, nSamples))
    for i in range(1,nOTUs+1):
        OTU_table_counts[i-1,:] = OTU_table_lines[i].split('\t')[1:]
    
    # Normalize counts by sample (column)
    for i in range(nSamples):
        OTU_table_counts[:,i] = np.divide(OTU_table_counts[:,i], np.sum(OTU_table_counts[:,i]))
    
    # Get sample contents for each taxa of the chosen level
    taxa_dict = {}
    for key in taxa_indices:
        taxa_dict[key] = {}
        indices = taxa_indices[key]
        for i in range(nSamples):
            taxa_dict[key][sample_labels[i]] = np.sum(OTU_table_counts[indices,i])

    return taxa_dict


def extract_available_taxa(OTU_table):
    # Extracts from a latin name OTU table the full list of taxa, from species to domain.  This can be used as a look-up for iterative calling of collapse_taxonomic_contents() and Rank Sum tests.
    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]
    
    # Taxonomic levels
    levels = ['domain',
              'phylum',
              'class',
              'order',
              'family',
              'genus',
              'species']

    # Create lists of each taxon from each level, its level and its corresponding parent
    taxa = []
    taxa_levels = []
    taxa_parents = []

    for taxonomic_level in levels:

        # Collapse to the right level        
        if(taxonomic_level == "domain"):
            OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "phylum"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[0][3:] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "class"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[1][3:] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "order"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[2][3:] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "family"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[3][3:] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "genus"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[4][3:] for OTU_ID in OTU_IDs]
        if(taxonomic_level == "species"):
            OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]
            OTU_taxa_parents = [OTU_ID.split(';')[5][3:] for OTU_ID in OTU_IDs]
        
        # Append to lists
        unique_set, unique_inds = np.unique(OTU_taxa, return_index=True)
        
        for i in range(len(unique_set)):
            taxa.append(OTU_taxa[unique_inds[i]])
            taxa_levels.append(taxonomic_level)
            if(taxonomic_level != "domain"):
                taxa_parents.append(OTU_taxa_parents[unique_inds[i]])
            else:
                taxa_parents.append('NA')

    return taxa, taxa_levels, taxa_parents


def get_taxonomic_level_name(latin_name_string, taxonomic_level):
    # Returns the chosen taxonomic level from an RDP type latin name string.  e.g. k__Bacteria;p__Firmicutes;..., 'phylum', returns 'Firmicutues'.
    latstr = latin_name_string
    if(taxonomic_level == "domain"):
        return latstr.split(';')[0].split('__')[-1]
    if(taxonomic_level == "phylum"):
        return latstr.split(';')[1].split('__')[-1]
    if(taxonomic_level == "class"):
        return latstr.split(';')[2].split('__')[-1]
    if(taxonomic_level == "order"):
        return latstr.split(';')[3].split('__')[-1]
    if(taxonomic_level == "family"):
        return latstr.split(';')[4].split('__')[-1]
    if(taxonomic_level == "genus"):
        return latstr.split(';')[5].split('__')[-1]
    if(taxonomic_level == "species"):
        return latstr.split(';')[6].split('__')[-1]


def extract_abundance_of_unassigned_OTUs(OTU_table, taxon_name):
    # Extracts the abundances of all OTUs that have a given taxonomy, e.g. 'k__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__'.

    # Load list of samples in the order they appear in the OTU table
    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        samples = OTU_table_lines[0].split('\t')[1:]
        samples[-1] = samples[-1].rstrip('\n')
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]

    # Get indices of each OTU corresponding to the desired taxon
    taxon_indices = []
    for i in range(len(OTU_taxa)):
        if OTU_taxa[i] == taxon_name:
            taxon_indices.append(i)

    if taxon_indices == []:
        return []

    # Read in OTU abundances
    nOTUs = len(OTU_IDs)
    nSamples = len(samples)

    # Load OTU table row major order
    OTU_table_counts = np.zeros((nOTUs, nSamples))
    for i in range(1,nOTUs+1):
        OTU_table_counts[i-1,:] = OTU_table_lines[i].split('\t')[1:]

    # Normalize counts by sample (column)
    for i in range(nSamples):
        OTU_table_counts[:,i] = np.divide(OTU_table_counts[:,i], np.sum(OTU_table_counts[:,i]))

    # Get sample contents for each taxa of the chosen level
    taxon_abundances = []
    for i in range(nSamples):
        taxon_abundances.append(np.sum(OTU_table_counts[taxon_indices, i]))

    return taxon_abundances


def extract_abundance_of_taxon(OTU_table, taxon_name, taxonomic_level, parent_node=0):
    # Returns a vector of abundances of the given taxon, where indices correspond to samples in the OTU table.  Takes as input the OTU table, the name of the taxon (e.g. 'Clostridium'), the taxonomic level (e.g. 'genus'), and the parent node (e.g. 'Clostridiaceae').

    # Load list of samples in the order they appear in the OTU table
    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        samples = OTU_table_lines[0].split('\t')[1:]
        samples[-1] = samples[-1].rstrip('\n')

    if taxonomic_level != "domain":
        taxon_dict = collapse_taxonomic_contents(OTU_table, taxonomic_level, parent_node)
        # Find the input taxon in the dict
        abundances = []
        for taxon in taxon_dict:
            single_name = get_taxonomic_level_name(taxon, taxonomic_level)
            if single_name == taxon_name:
                abundances = taxon_dict[taxon]
                break
        # populate vector in the order in which the samples come
        if len(abundances) > 0:
            abundance_vector = []
            for i in range(len(samples)):
                abundance_vector.append(abundances[samples[i]])
            return abundance_vector
        else:
            return 0
    else:
        print "Domain level not supported!"
        return 0


def collapse_taxonomic_contents_for_json(OTU_table, taxonomic_level, parent_node):
    # Takes as input an OTU table with latin names as OTU IDs,
    # the taxonomic level to which all taxonomies should be collapsed, and the parent taxon.
    # Accepts: domain, phylum, class, order, family, genus, species, subspecies as levels.
    # e.g.  collapse_taxonomic_contents(OTU_table, phylum, Bacteria)

    taxa_dict = collapse_taxonomic_contents(OTU_table, taxonomic_level, parent_node)
    sample_labels = []
    taxon_labels = []
    taxon_names = []

    # Get taxon names and sample IDs
    for taxon in taxa_dict:
        taxon_labels.append(taxon)
        taxon_names.append(taxon.split('__')[-1])
    try:
        for sampleID in taxa_dict[taxon_labels[-1]]:
            sample_labels.append(sampleID)
    except:
        print "No taxa found in samples."
    sorted_samples = sorted(sample_labels)
    sorted_sample_inds = [i[0] for i in sorted(enumerate(sample_labels), key=lambda x:x[1])]
    sorted_taxa = sorted(taxon_names)
    sorted_taxon_inds = [i[0] for i in sorted(enumerate(taxon_names), key=lambda x:x[1])]

    # Sort frequencies of each taxon for each sample
    frequencies = {}
    for sampleID in sorted_samples:
        frequencies[sampleID] = [taxa_dict[taxon_labels[index]][sampleID] for index in sorted_taxon_inds]
    frequencies['Taxa'] = sorted_taxa
    return frequencies


def extract_taxon_rankings(OTU_table, taxonomic_level, parent_node=0):
    # Takes as input an OTU table with latin names.
    # Returns a dict with sampleIDs as keys: dict[sampleID]['names'] = [clade1, clade2, clade3], dict[sampleID]['abundances'] = [abundance1, abundance2, abundance3]  
    clade_dict = collapse_taxonomic_contents(OTU_table, taxonomic_level, parent_node)
    output_dict = {}
    clades = clade_dict.keys()
    sampleIDs = clade_dict[clades[0]].keys()
    # For each sample, populate a vector with clade abundances
    for sampleID in sampleIDs:
        output_dict[sampleID] = {}
        output_dict[sampleID]['names'] = []
        output_dict[sampleID]['abundances'] = []
        clade_vector = []
        for clade in clades:
            clade_vector.append(clade_dict[clade][sampleID])
        # Sort vector
        sorted_inds = np.argsort(clade_vector)[::-1]
        for ind in sorted_inds:
            output_dict[sampleID]['names'].append(clades[ind])
            output_dict[sampleID]['abundances'].append(clade_vector[ind])
    return output_dict
    

