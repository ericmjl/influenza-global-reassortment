from collections import Counter

import numpy as np
import custom_funcs as cf

"""
This Python module contains the accessory functions necessary to compute the data for figure plotting. 

Do not delete this module, as it contains crucial business logic for plotting the data. Only functions that are reused
across notebooks are saved here.
"""

# This set of functions are used to compute the data needed to plot the known transmission patterns.
def same_host_descent(G, host_species):
    """
    Returns the nodes that are involved in same-host transmission.
    """
    nodes = set()  # nodes that are involved in human-human transmission.
    for sc, sk, d in G.edges(data=True):
        sc_host = G.node[sc]['host_species']
        sk_host = G.node[sk]['host_species']
        
        sc_subtype = G.node[sc]['subtype']
        sk_subtype = G.node[sk]['subtype']
        
        subtype_exclusions = ['Mixed', 'mixed']
        
        not_mixed = (sc_subtype not in subtype_exclusions) and (sk_subtype not in subtype_exclusions)
        is_clonal = d['edge_type'] == 'full_complement'
        
        if sc_host == host_species and sk_host == host_species and not_mixed and is_clonal:
            nodes.add(sc)
            nodes.add(sk)

    return nodes

def subtype_counts(node_set, G, log=False):
    """
    Returns a dictionary of subtypes and the number of nodes that are present.
    
    Parameters:
    - node_set: a list of nodes
    - G: a networkX graph from which the nodes come from. 
    """
    subtypes = Counter()
    for n in node_set:
        subtype = G.node[n]['subtype']
        subtypes[subtype] += 1

    if log:
        for k, v in subtypes.items():
            subtypes[k] = np.log10(v)
            
    return subtypes


### Used for making proportion reassortant under null distribution. ###
# Computes a null distribution proportion for host species

def null_distribution_proportion_reassortant(G, equally=False):
    
    excluded_hosts = ['Aquatic Bird', 'Avian', 'Bird', 'Duck', 'Environment', 'Mallard-Black Duck Hybrid', 'Sea Mammal', 
                  'Unknown', 'Waterfowl']
    
    G_shuffled = cf.shuffle_node_attribute_label(G, 'host_species', equally)
    props = cf.edge_proportion_reassortant(G_shuffled, 'host_species', exclusions=excluded_hosts)
    return props