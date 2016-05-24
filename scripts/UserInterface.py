from ete2 import Tree
import json
import Taxonomy
import os

def taxonomy_parser(node_label):
# Taxonomy parser for 'fill_mod.newick' tree
    taxa = ['domain','phylum','class','order','family','genus','species','subspecies','strain']
    labels = node_label.split('_')
    taxonomic_properties = {}
    labels = [x for x in labels if x != "diff"]
    if labels[0] != "Root":
        for i in range(len(labels)):
            taxonomic_properties[taxa[i]] = labels[i]
        taxonomic_properties['Level'] = taxa[i]
        taxonomic_properties['Children'] = taxa[i+1]
    else:
        taxonomic_properties['Level'] = 'Root'
        taxonomic_properties['Children'] = taxa[0]
    return taxonomic_properties

def write_json_files(OTU_table, target_directory):
    # Write JSON files containing children node abundances from a given OTU table to a target directory.
    # Currently does this for each internal node in the 'fill_mod.newick' tree.
    tree_file = '/home/ubuntu/templates/fill_mod.newick'
    OTU_table_labeled = OTU_table + '.taxonomies'
    Taxonomy.convert_GGIDs_to_latin_names(OTU_table, OTU_table_labeled, '/home/ubuntu/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt')

    # Load tree
    tree = Tree(tree_file, format=1)
    leaves = []
    internals = []
    all_nodes = tree.get_descendants("preorder")

    # Sort into leaves and internal nodes
    for node in all_nodes:
        if(node.is_leaf()):
            leaves.append(node)
        else:
            internals.append(node)

    # Write each node's JSON
    for node in internals:
        node_taxonomies = taxonomy_parser(node.name)
        node_name = node_taxonomies[node_taxonomies['Level']]
        children_taxa = node_taxonomies['Children']
        json_dict = Taxonomy.collapse_taxonomic_contents_for_json(OTU_table_labeled, children_taxa, node_name)
        with open(os.path.join(target_directory, node_name + '.json'), 'w') as outfile:
            if len(json_dict) > 0:
                json.dump(json_dict, outfile, sort_keys=True)

