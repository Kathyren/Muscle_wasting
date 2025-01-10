import pandas as pd
from networkx import is_bipartite
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cytoscape.enum_network_sources import NetworkSource

import Constants
import networkx as nx
from network import network_processing
from network.main_pipeline import get_cytoscape_network
from network.network_processing import create_graph, draw_graph, load_graph, \
    create_graph_from_dictionaries, get_mirna_mrna_relationships, remove_nodes_low_centrality, convert_to_json, \
    save_graph, get_nodes_names, add_mirna_relationships, set_positions, get_mirna_tissue_edges, \
    add_tissue_relationship, add_organ_system_relationship, get_tissue_system_edges, extract_genes_from_pathways, \
    add_pathway_to_node, add_pathways_to_nodes, mark_TF_nodes_from_file, add_dds_to_node, add_DDS_data, random_walk
from cytoscape import read_cytoscape_json, format_cytoscape_json, create_cytoscape_node, create_cytoscape_edge, protein_name

import pytest


def test_use_example():
    """
    This function will use the example to create a network. It is not a proper test, is to 
    remmember me how to use it. 
    TBD 
    This will be what is expected to be in the manual.
    :return:
    """

    cytoscape_json = read_cytoscape_json(cytoscape_file='test.cyjs')
    # graph = load_graph()
    c_nodes, c_edges, c_relationships = format_cytoscape_json(cytoscape_json)
    graph = create_graph_from_dictionaries(nodes=c_nodes, relationship=c_relationships, edges=c_edges)
    # graph = create_my_graph(get_random_relationships)
    save_graph(graph, "small_graph.pkl")
    anc = nx.all_pairs_node_connectivity(graph)

    print(nx.info(graph))

    # draw_graph(graph)
