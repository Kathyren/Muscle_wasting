### This file is the mail pipeline for my priject
import cytoscape as ct
import networkX as nx
import os.path
from os import path
import py4cytoscape as py4


def get_cytoscape_network(file_name):
    cytoscape_json = ct.read_cytoscape_json(cytoscape_file=file_name)
    nodes, edges_metadata, relationship = ct.format_cytoscape_json(cytoscape_json=cytoscape_json)
    the_network = nx.create_graph_from_dictionaries(nodes=nodes, edges=edges_metadata, relationship=relationship)
    return the_network


def main(file, name):
    if not os.path.exists(f"graph0_{name}"):
        the_network = get_cytoscape_network(file)
        nx.save_graph(the_network, f"graph0_{name}")
    else:
        the_network = nx.load_graph(f"graph0_{name}")
    network = nx.remove_nodes_low_centrality(graph=the_network, cutoff=0.75)
    save_as_cjsn(network, f'graph1_{name}')


def save_as_cjsn(network, name):
    """
    This function will take a network and save it for cytoscape use
    :return:
    """
    json_network = nx.convert_to_json(network)
    ct.save_cytoscape_json(json_network, cytoscape_file_name=name)


def open_cytoscape(file_name):
    py4.open_session(file_name)


if __name__ == '__main__':
    file = ""
    main(file, "test")
    pass
