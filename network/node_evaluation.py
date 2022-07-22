import networkx as nx



def calculate_all_nodes_centralities(G):
    """
    This function will calculate the
    -In degree centrality
    -out degree centrality
    -betweenness centrality
    -Load centralitu
    -Eigenvector centrality
    -closeness centrality
    Depending on the use (they are commented now)

    :param G: A networkX graph
    :return: a list of list of the centrality for each node.
    """
    centralities = []
    # dc = list(nx.degree_centrality(graph).values())
    # idc = nx.in_degree_centrality(G)
    # odc = nx.out_degree_centrality(G)
    # bc = list(nx.betweenness_centrality(graph).values())
    # lc = list(nx.load_centrality(graph).values())
    #    ec = list(nx.eigenvector_centrality(G).values())
    cc_t = nx.closeness_centrality(G)
    cc = list(cc_t.values())
    centralities.append(cc)
    return centralities