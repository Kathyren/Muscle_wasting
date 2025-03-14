import sys
sys.path.append('../')
sys.path.append('../network')
from mirna_scoring.walking_network import *
import networkx as nx
import pandas as pd

def test_traverse_and_update():
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/sub_network_nodes.pkl'
    graph = nx.read_gpickle(graph_pkl)
    strat_node = 'hsa-miR-21-5p'
    a, b, c = traverse_and_update(graph, strat_node)

    print(a)
    print(b)
    print(c)
    d = {'DNM1L': [-1],
         'APAF1': [-1],
         'LRRFIP1': [-1],
         'IL1B': [1, -1],
         'CAMK2N1': [-1],
         'E2F3': [-1],
         'GAPDH': [-1],
         'MAP3K5': [-1],
         'MYC': [-1],
         'PPIF': [-1, -1],
         'SMAD7': [-1, -1],
         'TPM1': [-1, 1],
         'VEGFA': [-1, -1, -1],
         }
def test_start_mir_path():
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
    graph = nx.read_gpickle(graph_pkl)
    strat_node = 'hsa-miR-122-5p'
    start_mir_path(graph, strat_node)

    nodes = graph.nodes
    #import network.network_processing as np
    #np.save_graph(graph, f"/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl")

    ALDOA = nodes['ALDOA']['data']['influence'][strat_node]
    PKM = nodes['PKM']['data']['influence'][strat_node]
    assert ALDOA==[-1]
    assert PKM==[-1]
    assert 'influence' not in nodes['NT5E']['data'], f"No mirna influences NT5E"

def test_evaluate_pathway_influence():
    influence_data = pd.DataFrame(columns=['ALDOA', 'PKM'], data=[[[-1], [-1]]], index=['hsa-miR-122-5p']).T
    mir = 'hsa-miR-122-5p'
    influence_data = [{'m_l': 0, 'm_s': 0, 'mo': 0, 'pathways': ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', 'KEGG_PYRUVATE_METABOLISM'],'ym': 0, 'yo': 0},
                      {'m_l': 0, 'm_s': 0, 'mo': 0, 'pathways': ['GOBP_SMALL_MOLECULE_METABOLIC_PROCESS', 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS'], 'ym': 0, 'yo': 0}]
    pass

def test_evaluate_pathway_influence():
    pass
def test_register_path():
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/sub_network_nodes.pkl'
    graph = nx.read_gpickle(graph_pkl)
    strat_node = 'hsa-miR-145-5p'
    node = graph.nodes[strat_node]
    node['data']['influence'] = {strat_node: [1]}
    p = register_path(graph, node, strat_node, visited_edges=[])
    print(p)

def test_evaluate_path_influence():

    pass
def initialize_invariant(graph):
    for node in graph.nodes:
        assert graph.nodes[node]['data']['influence'] == [1], "Initialization invariant failed"

def traversal_invariant(node, edge, weight, visited_edges):
    assert edge not in visited_edges, "Traversal invariant failed: edge visited more than once"
    assert weight == edge['data']['weight'], "Traversal invariant failed: weight mismatch"

def update_invariant(node, influence_value):
    assert node['data']['influence'][-1] == influence_value, "Update invariant failed"

def termination_invariant(graph, visited_nodes, initial_nodes):
    for node in graph.nodes:
        if node in visited_nodes or node in initial_nodes:
            assert 'influence' in graph.nodes[node]['data'], "Termination invariant failed: node influence not updated"
