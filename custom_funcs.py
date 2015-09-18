from random import shuffle

# For cleaning host species names.
def clean_host_species_names(G):
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

def shuffle_node_attribute_label(G, attribute):
    # Some initial checks
    for n, d in G.nodes(data=True):
        try:
            assert attribute in d.keys()
        except AssertionError:
            print('Attribute {0} not in node {1} metadata'.format(attribute, n))

    attrs = [d[attribute] for n, d in G.nodes(data=True)]
    shuffle(attrs)

    G_shuffled = G.copy()

    for i, (n, d) in enumerate(G_shuffled.nodes(data=True)):
        G_shuffled.node[n][attribute] = attrs[i]

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

            if sk_attr not in exclusions and sc_attr not in exclusions:
                if sk_attr == sc_attr:
                    n_same += 1
                else:
                    n_diff += 1

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
