import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import warnings
import numpy as np

from database_analysis import sql_operations as sql


def create_graph(mirnas, genes, relationsip):
    """
    This function will receive a list of genes and mirnas and create the network of them
    :param mirnas: List of string with the mirna names
    """
    G = nx.Graph()
    G.add_nodes_from(mirnas, size=10, type="mirna")
    G.add_nodes_from(genes, type="gene")
    G.add_edges_from(relationsip)
    is_bipartite(G)
    return G


def draw_graph(G):
    color_map = []

    for node, data in G.nodes(data=True):
        if data['type'] == 'mirna':
            color_map.append(0.25)  # blue color
        elif data['type'] == 'gene':
            color_map.append(0.7)  # yellow color
    nx.draw(G, node_color=color_map, with_labels=False)
    # nx.draw(G)


def is_bipartite(G):
    return True
    for edge, data in G.edges(data=True):
        count = 0
        for node in edge:
            if '-miR-' in node:
                count = count + 1
        assert count == 1, 'There is a mirna_mirna interaction or a gene- gene interaction'
    return True

def filter_by_degree(G):
    graph = G.copy()
    centrality = nx.degree(G)
    c = list(centrality)
    _, values = zip(*c)
    average = np.mean(values)
    sandar_dev = np.std(values)
    x = np.quantile(values, 0.75)
    delete_nodes = []
    for node in centrality:
        if node[1] < x:
            delete_nodes.append(node[0])
    graph.remove_nodes_from(delete_nodes)
    save_graph(graph, "high_degree.pkl")
    return  graph
def calculate_centralities(G):
    graph = load_graph("high_degree.pkl")

    #dc = list(nx.degree_centrality(graph).values())
    # idc = nx.in_degree_centrality(G)
    # odc = nx.out_degree_centrality(G)
    #bc = list(nx.betweenness_centrality(graph).values())
    #lc = list(nx.load_centrality(graph).values())
    #    ec = list(nx.eigenvector_centrality(G).values())
    cc_t = nx.closeness_centrality(graph)
    cc = list(cc_t.values())

    x = np.quantile(cc, 0.75)
    delete_nodes = []
    for node in cc_t.items():
        if node[1] < x:
            delete_nodes.append(node[0])
    graph.remove_nodes_from(delete_nodes)
    draw_graph(graph)
    pass

def calculate_data(dc, bc, lc, cc):
    centralities = [dc,bc,lc, cc]
    av = []
    std = []
    max = []
    min = []
    for c in centralities:
        av.append(np.mean(c))
        std.append(np.std(c))
        max.append((np.max(c)))
        min.append(np.min(c))
    print(av, std, max, min)
def save_graph(G, name):
    nx.write_gpickle(G, name)


def load_graph(name):
    G = nx.read_gpickle(name)
    return G


get_random_relationships = "Select Distinct mrna, mirna_mature from binding limit 10;"
get_muscle_relationships = 'Select Distinct mrna, mirna_mature from binding where  ' \
                           'mirna_mature like "hsa-%" and mrna in (select gene from gene_bank) ;'


def create_my_graph(query):
    relationships = sql.get_query(query=query)
    genes = list(relationships['mrna'])
    mirnas = list(relationships['mirna_mature'])
    relationship = list(zip(genes, mirnas))
    graph = create_graph(mirnas=mirnas, genes=genes, relationsip=relationship)
    return graph


if __name__ == '__main__':
    # graph = load_graph()
    graph = create_my_graph(get_random_relationships)
    save_graph(graph, "small_graph.pkl")

    # draw_graph(graph)
