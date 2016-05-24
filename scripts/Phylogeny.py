"""

OVERVIEW: 

Module to compute phylogenetic quantities.

"""

import os
import Formatting as frmt
from ete2 import Tree
import math

def ComputeTreeAbundances(OTU_table_file, tree_file):
        # Computes L-R type feature from a tree and OTU table in classic dense format.
        # NOTE: does not support GreenGenes tree yet (e.g. '97_otus.tree') due to internal node nomenclature being non-standard.  Needs quick fix.

        # Load abundances
        OTUs, samples, OTU_table = frmt.load_OTU_table_classic(OTU_table_file)
        abundances = {}
        for i in range(len(OTUs)):
            abundances[OTUs[i]] = {}
            for j in range(len(samples)):
                abundances[OTUs[i]][samples[j]] = OTU_table[i,j]
		
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

        # Add root node to internal node list
        ancestors = [node for node in leaves[0].iter_ancestors()]
        root = ancestors[len(ancestors)-1]
        internals.append(root)

        # Get the parent node of each node
        parent_map = {}
        for node in internals:
            # For each internal node, fill its children's parent entry
            children = node.get_children()
            parent_map[children[0]] = node
            parent_map[children[1]] = node


        # Each internal node has dictionary entries of type nodes[node][L] and nodes[node][R].
        # For each leaf node, add it to each successive parent node working up the tree, either as left or right
        nodes = {}
        for leaf in leaves:
            newleaf = leaf
            ancestor_iterator = leaf.iter_ancestors()
            for ancestor in ancestor_iterator:
                if ancestor.name not in nodes:
                        nodes[ancestor.name] = {}
                        nodes[ancestor.name]['L'] = []
                        nodes[ancestor.name]['R'] = []
                children = ancestor.get_children()
                if (newleaf == children[0]):
                        nodes[ancestor.name]['L'].append(leaf.name)
                elif (newleaf == children[1]):
                        nodes[ancestor.name]['R'].append(leaf.name)
                newleaf = ancestor


        # Compute L-R for each internal node
        LminusR = {}

        for node in nodes:
            LminusR[node] = {}
            for sampleID in samples:
                left_nodes = nodes[node]['L']
                left_abundances = sum([abundances[leftnode][sampleID] for leftnode in left_nodes])
                right_nodes = nodes[node]['R']
                right_abundances = sum([abundances[rightnode][sampleID] for rightnode in right_nodes])
                try:
                        left_abundances = math.log(left_abundances)
                except:
                        left_abundances = -16
                try:
                        right_abundances = math.log(right_abundances)
                except:
                        right_abundances = -16
                LminusR[node][sampleID] = left_abundances - right_abundances

        return LminusR

