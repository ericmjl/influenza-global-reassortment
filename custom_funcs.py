from random import shuffle, choice
import pandas as pd
from itertools import combinations

def clean_host_species_names(G):
    """
    Cleans host species names according to rules specified.

    Parameters:
    ===========
    - G: (nx.DiGraph) The network of viral isolates.

    Returns:
    ========
    - G: (nx.DiGraph) The network of viral isolates, with host species names 
         cleaned.
    """
    for n, d in G.nodes(data=True):
        host_species = d['host_species']
        if '/' in host_species:
            host_species = host_species.split('/')[0]
            G.node[n]['host_species'] = host_species

        if 'Mallard Duck' in host_species:
            G.node[n]['host_species'] = 'Mallard'
        if 'Eurasian' in host_species:
            G.node[n]['host_species'] = host_species.split(' ')[1]
        if 'Pintail' in host_species:
            G.node[n]['host_species'] = 'Pintail'


    return G

def impute_weights(G):
    """
    This function exists to impute weights. 

    If a node is a clonal node, but has multiple sources, then it should be 
    
        1 / n_sources 
    
    weight per edge.

    If a node is reassortant, but has more than 1 possible source pair, i.e.
    source 1 and 2, or source 1 and 3, then the double-counted edge should be
    weighted heavier than the single-counted edge. Mutadis mutandis for
    triple, quadruple etc.
    """
    ### SOME CHECKS BEFORE RUNNING FUNCTION ###
    try:
        for n, d in G.nodes(data=True):
            assert 'reassortant' in d.keys()
    except AssertionError:
        print('Impute reassortant status before calling this function.')

    ### ACTUAL FUNCTION BEGINS BELOW ###
    for n, d in G.nodes(data=True):
        
        # Initialize weight counts
        for sc, sk, _ in G.in_edges(n, data=True):
            G.edge[sc][sk]['weight_ct'] = 0
        
        if d['reassortant']:
            for e1, e2 in combinations(G.in_edges(n, data=True), 2):
                sc1, sk1, d1 = e1
                sc2, sk2, d2 = e2
                
                # Grab out the segments.
                segs1 = d1['segments'].keys()
                segs2 = d2['segments'].keys()

                if set(segs1).union(segs2) == set(range(1,9)):
                    G.edge[sc1][sk1]['weight_ct'] += 1
                    G.edge[sc2][sk2]['weight_ct'] += 1

        if not d['reassortant']:
            for sc, sk, d in G.in_edges(n, data=True):
                G.edge[sc][sk]['weight_ct'] += 1

        weight_ct_total = sum([d['weight_ct'] for sc, sk, d in G.in_edges(n, data=True)])

        for sc, sk, d in G.in_edges(n, data=True):
            weight_ct = d['weight_ct']
            if weight_ct_total > 0:
                G.edge[sc][sk]['weight'] = weight_ct / weight_ct_total
            else:
                G.edge[sc][sk]['weight'] = 0

    return G

def remove_zero_weighted_edges(G):
    for sc, sk, d in G.edges(data=True):
        if d['weight'] == 0:
            G.remove_edge(sc, sk)
    return G

def impute_reassortant_status(G):
    for n, d in G.nodes(data=True):
        in_edges = G.in_edges(n, data=True)
        for sc, sk, d in in_edges:
            if d['edge_type'] == 'reassortant':
                G.node[n]['reassortant'] = True
                break
            else:
                G.node[n]['reassortant'] = False
                break
        if len(in_edges) == 0:
            G.node[n]['reassortant'] = False
    return G

def shuffle_node_attribute_label(G, attribute, equally=False):
    # Some initial checks
    for n, d in G.nodes(data=True):
        try:
            assert attribute in d.keys()
        except AssertionError:
            print('Attribute {0} not in node {1} metadata'.format(attribute, n))

    attrs = [d[attribute] for n, d in G.nodes(data=True)]
    if equally == False:
        shuffle(attrs)

    elif equally == True:
        attrs_set = list(set(attrs))

    G_shuffled = G.copy()

    for i, (n, d) in enumerate(G_shuffled.nodes(data=True)):
        if equally == False:
            G_shuffled.node[n][attribute] = attrs[i]
        elif equally == True:
            G_shuffled.node[n][attribute] = choice(attrs_set)

    return G_shuffled

def count_edges(G, attr, exclusions=[]):
    """
    Returns a dictionary of edge counts for full complement vs. reassortant in 
    same hosts and different hosts, such that all counts are normalized to 1 
    for each node.
    """

    edge_counts = dict()
    edge_counts['full_complement'] = dict()
    edge_counts['full_complement']['same_attr'] = 0
    edge_counts['full_complement']['diff_attr'] = 0
    edge_counts['reassortant'] = dict()
    edge_counts['reassortant']['same_attr'] = 0
    edge_counts['reassortant']['diff_attr'] = 0


    for n, node_d in G.nodes(data=True):
        in_edges = G.in_edges(n)
        n_same = 0
        n_diff = 0
        is_reassortant = node_d['reassortant']
        sk_attr = G.node[n][attr]

        for sc, _, edge_d in G.in_edges(n, data=True):
            sc_attr = G.node[sc][attr]
            weight = edge_d['weight']

            if sk_attr not in exclusions and sc_attr not in exclusions:
                if sk_attr == sc_attr:
                    n_same += weight
                else:
                    n_diff += weight

        # Calculate the relative proportion of each, else set proportions to 0.
        if n_same + n_diff > 0:
            p_same = n_same / (n_same + n_diff)
            p_diff = n_diff / (n_same + n_diff)

        if n_same + n_diff == 0:
            p_same = 0
            p_diff = 0

        if is_reassortant:
            edge_counts['reassortant']['same_attr'] += p_same
            edge_counts['reassortant']['diff_attr'] += p_diff
        else:
            edge_counts['full_complement']['same_attr'] += p_same
            edge_counts['full_complement']['diff_attr'] += p_diff

    return edge_counts

def edge_proportion_reassortant(G, attr, exclusions=[]):
    edge_counts = count_edges(G, attr, exclusions)

    same_attr = edge_counts['reassortant']['same_attr'] / (edge_counts['reassortant']['same_attr'] + edge_counts['full_complement']['same_attr'])

    diff_attr = edge_counts['reassortant']['diff_attr'] / (edge_counts['reassortant']['diff_attr'] + edge_counts['full_complement']['diff_attr'])

    return {'same_attr':same_attr, 'diff_attr':diff_attr}

def impute_host_group_name(G):
    hgs = pd.read_csv('host_groups.csv')
    hgs_dict = dict()
    for r, d in hgs.iterrows():
        hgs_dict[(d['Country'], d['Species'])] = d['Habitat/setting']
        
    for n, d in G.nodes(data=True):
        host = d['host_species']
        country = d['country']
        G.node[n]['host_group'] = hgs_dict[(country, host)]

    return G