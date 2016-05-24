"""

OVERVIEW: 

Python module for microbiome statistical analysis tools.

"""

import os, sys
import numpy as np
import warnings
import scipy.stats
import Taxonomy
    
def jsd(x,y): 
    #Jensen-shannon divergence
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log2(2*x/(x+y))
    d2 = y*np.log2(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)    
    return d


def taxon_ranksumtests(OTU_table_1, OTU_table_2, outputfile, taxonomies = '/home/ubuntu/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt'):
    # Performs Wilcoxon Rank Sum tests between two OTU tables with the same samples, for every possible clade and OTU based on GreenGenes references    
    OTU_table_1_latin = OTU_table_1 + '.latin'
    OTU_table_2_latin = OTU_table_2 + '.latin'
    Taxonomy.convert_GGIDs_to_latin_names(OTU_table_1, OTU_table_1_latin, taxonomies)
    Taxonomy.convert_GGIDs_to_latin_names(OTU_table_2, OTU_table_2_latin, taxonomies)

    # Get list of samples
    with open(OTU_table_1, 'r') as fid:
        all_lines = fid.readlines()
        samples = all_lines[0].split('\t')[1:]
        samples[-1] = samples[-1].rstrip('\n')

    # Load OTU abundances
    with open(OTU_table_1, 'r') as fid:
        disease_lines = fid.readlines()[2:]
    with open(OTU_table_2, 'r') as fid:
        healthy_lines = fid.readlines()[2:]

    OTU_IDs = [line.split()[0] for line in disease_lines]
    OTU_table_1_dict = {}
    OTU_table_2_dict = {}
    for line in disease_lines:
        split_line = line.split()
        x = np.array(split_line[1:], dtype='|S4')
        vals = x.astype(np.float)
        OTU_table_1_dict[split_line[0]] = vals
    for line in healthy_lines:
        split_line = line.split()
        x = np.array(split_line[1:], dtype='|S4')
        vals = x.astype(np.float)
        OTU_table_2_dict[split_line[0]] = vals

    # Get list of taxa and loop through this list to get abundances; add these abundance vectors to the dicts above
    taxa, taxa_levels, taxa_parents = Taxonomy.extract_available_taxa(OTU_table_1_latin)

    ntaxa = len(taxa)
    for i in range(ntaxa):
        level = taxa_levels[i]
        parent = taxa_parents[i]
        taxonomic_level_name = Taxonomy.get_taxonomic_level_name(taxa[i], level)
        if parent != "NA" and len(taxonomic_level_name) > 0:
            abundances_table_1 = Taxonomy.extract_abundance_of_taxon(OTU_table_1_latin, taxonomic_level_name, taxa_levels[i], taxa_parents[i])
            abundances_table_2 = Taxonomy.extract_abundance_of_taxon(OTU_table_2_latin, taxonomic_level_name, taxa_levels[i], taxa_parents[i])
            if abundances_table_1 !=  0 and abundances_table_2 != 0:
                OTU_table_1_dict[taxa[i]] = abundances_table_1
                OTU_table_2_dict[taxa[i]] = abundances_table_2

    # Wilcoxon Rank Sum test
    p_vals = []
    z_stats = []
    taxa_names = []
    for key in OTU_table_1_dict:
        if OTU_table_1_dict[key] is not None and OTU_table_2_dict[key] is not None:
            [zstat, pval] = scipy.stats.ranksums(OTU_table_1_dict[key], OTU_table_2_dict[key])
            taxa_names.append(key)
            p_vals.append(pval)
            z_stats.append(zstat)


    # Use Bonferroni correction for multiple testing
    alpha = 0.05/len(p_vals)

    # Find p-values less than significance cut-off
    significant_indices = []
    counter = 0
    for p_val in p_vals:
        if p_val < alpha:
            significant_indices.append(counter)
        counter += 1


    # Order by p-value and write to file
    sorted_inds = np.argsort(p_vals)

    # Write to file
    with open(outputfile, 'w') as fid:
        firstline = '\t'.join(['Wilcoxon P-value', 'Significant','Taxon', 'GGID'])
        fid.write(firstline + '\n')
        for index in sorted_inds:
            significant = 'False'
            if index in significant_indices:
                significant = 'True'
            if taxa_names[index][:3] != "k__":
                taxon_name = Taxonomy.GGID_lookup(taxa_names[index])
                GGID = taxa_names[index]
            else:
                taxon_name = taxa_names[index]
                GGID = 'N/A'
            line = '%s' % float('%.3g' % p_vals[index]) + '\t' + significant + '\t' + str(taxon_name) + '\t' + str(GGID)
            fid.write(line + '\n')

