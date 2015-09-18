def clean_host_species_names(G):
    G = G.copy()
    for n, d in G.nodes(data=True):
        