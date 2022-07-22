import networkx as nx
import numpy as np


def distance_to_target(network, source, target='muscle'):
    """

    :param network:
    :param source: str Node name (commonly the mirna)
    :param target: str Node name (by now, the muscle)
    :return: list with the route
    """
    paths = list(nx.shortest_simple_paths(network, source=source, target=target, weight='weight'))
    print(paths)
    shortest_paths = nx.shortest_path(network, source=source, target=target, weight='weight')
    return shortest_paths


def evaluate_nodes(network):
    mirna_nodes = [x for x, y in network.nodes(data=True) if y['type'] == 'mirna']
    ec = nx.eigenvector_centrality(network, weight='weight', tol=1.0e-3)
    for mirna in mirna_nodes:
        path = distance_to_target(network, source=mirna, target='muscle')
        weight = nx.path_weight(network, path=path, weight='weight')
        network.nodes[mirna]["shortest_path_len"] = len(path)
        network.nodes[mirna]["path_weight"] = weight
        centrality = ec[mirna]
        network.nodes[mirna]["centrality"] = centrality

def remove_nodes(network, threshold):
    mirna_nodes = [x for x, y in network.nodes(data=True) if y['type'] == 'mirna']
    total_scores = []
    alpha = 1
    beta = 1
    gama = 1
    for mirna in mirna_nodes:
        node = network.nodes[mirna]
        length = node["shortest_path_len"]
        weight = node["path_weight"]
        centrality = node["centrality"]
        total_score = beta*weight+gama*centrality - alpha*length
        total_scores.append(total_score)
    cutoff = np.quantile(total_scores, threshold)
    delete_nodes = []
    for idx, value in enumerate(total_scores):
        if value < cutoff:
            delete_nodes.append(mirna_nodes[idx])
    network.remove_nodes_from(delete_nodes)

def shortest_distance_max_value_to_target(network, source, target='muscle'):
    """

    :param network:
    :param source: str Node name (commonly the mirna)
    :param target: str Node name (by now, the muscle)
    :return: list with the route
    """
    paths = list(nx.shortest_simple_paths(network, source=source, target=target, weight='weight'))
    print(paths)
    path_lens = [len(path) for path in paths]
    path_weight = []
    for path in paths:
        weight = nx.path_weight(network, path=path, weight='weight')
        path_weight.append(weight)


def get_pareto(lengths, weights):
    """
    This function will retrieve the pareto set of the
    max year
    max cites
    :param plot: If we should plot the pareto
    :return:
    """
    length_weight = []
    for idx, length in enumerate(lengths):
        try:
            length_weight.append([length, -weights[idx]])
        except ValueError:
            length_weight.append([999, 0])
    npArray = np.array(length_weight)
    pareto = is_pareto_efficient(npArray)
    pareto_bool = pareto.tolist()
    return pareto_bool


def is_pareto_efficient(costs, return_mask=True):
    """
    Find the pareto-efficient points
    Obtained from https://github.com/QUVA-Lab/artemis/blob/peter/artemis/general/pareto_efficiency.py
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index < len(costs):
        nondominated_point_mask = np.any(costs < costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype=bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient


def calculate_all_nodes_closeness_centrality(graph):
    """
    This function will calculate the
    -closeness centrality

    :param graph: A networkX graph
    :return: a list of tuples node-centrality
    """
    cc_t = nx.closeness_centrality(graph)
    nodes = list(graph.nodes)
    cc = list(cc_t.values())
    centralities = zip(nodes, cc)
    return centralities


def calculate_all_nodes_betweenness_centrality(graph):
    """
    This function will calculate the
    -betweenness centrality

    :param graph: A networkX graph
    :return:  a list of tuples node-centrality
    """
    bc = list(nx.betweenness_centrality(graph).values())
    nodes = list(graph.nodes)
    cc = list(bc.values())
    centralities = zip(nodes, cc)
    return centralities


def calculate_all_nodes_eigenvector_centrality(graph, weights='weight'):
    """
    This function will calculate the
    -betweenness centrality

    :param weights: value in the network to consider as weight
    :param graph: A networkX graph
    :return: dictionary with the name of the node and the centralities
    """
    ec = nx.eigenvector_centrality(graph, weight=weights)
