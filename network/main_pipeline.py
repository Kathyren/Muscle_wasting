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


def add_mirnas_n_select(cytoscape_network, name, use_prefix=True, add_mirnas=True):
    if use_prefix:
        px1 = "graph0_"
        px2 = "graph1_"
    else:
        px1 = px2 = ""
    if not os.path.exists(f"{px1}{name}.pkl"):
        the_network = get_cytoscape_network(cytoscape_network)
        nx.save_graph(the_network, f"{px1}{name}.pkl")
    else:
        the_network = nx.load_graph(f"{px1}{name}.pkl")
    if add_mirnas:
        nx.add_mirna_relationships(the_network)

    network = nx.remove_nodes_low_centrality(graph=the_network, cutoff=0.5)
    nx.set_positions(network)
    nx.save_graph(network, f"{px2}{name}.pkl")
    save_as_cjsn(network, f'{px2}{name}.cyjs')


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
    file = "GSE38718.cyjs"
    magagnes2009 = "STRING_Magagnes2009"
    #add_mirnas_n_select(file, "dryrun_cardiovascular_w_mirnas")
    #add_mirnas_n_select(file, "GSE38718_w_mirnas")
    add_mirnas_n_select(magagnes2009+".cyjs", magagnes2009+"2")
    pass