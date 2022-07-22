import networkx as nx


def distance_to_target(network, source, target='muscle'):
    """

    :param source: str Node name (commonly the mirna)
    :param target: str Node name (by now, the muscle)
    :return: list with the route
    """

    shortest_paths = nx.shortest_path(network, source=source, target=target, weight='weight')
    return shortest_paths


def calculate_all_nodes_closeness_centrality(G):
    """
    This function will calculate the
    -closeness centrality

    :param G: A networkX graph
    :return: a list of tuples node-centrality
    """
    cc_t = nx.closeness_centrality(G)
    nodes = list(G.nodes)
    cc = list(cc_t.values())
    centralities = zip(nodes, cc)
    return centralities

def calculate_all_nodes_betweenness_centrality(G):
    """
    This function will calculate the
    -betweenness centrality

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
