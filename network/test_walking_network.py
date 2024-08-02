from walking_network import *
import networkx as nx


def test_traverse_and_update():
    graph_pkl = 'sub_network_nodes.pkl'
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
    graph_pkl = 'sub_network_nodes.pkl'
    graph = nx.read_gpickle(graph_pkl)
    strat_node = 'hsa-miR-21-5p'
    start_mir_path(graph, strat_node)

    nodes = graph.nodes
    PPIF = nodes['PPIF']['data']['influence'][strat_node]
    SMAD7 = nodes['SMAD7']['data']['influence'][strat_node]
    VEGFA = nodes['VEGFA']['data']['influence'][strat_node]
    TPM1 = nodes['TPM1']['data']['influence'][strat_node]

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

    assert VEGFA == d['VEGFA']
    assert TPM1 == d['TPM1']
    assert SMAD7 == d['SMAD7']
    assert PPIF == d['PPIF']
def test_register_path():
    graph_pkl = 'sub_network_nodes.pkl'
    graph = nx.read_gpickle(graph_pkl)
    strat_node = 'hsa-miR-145-5p'
    node = graph.nodes[strat_node]
    node['data']['influence'] = {strat_node: [1]}
    p = register_path(graph, node, strat_node, visited_edges=[])
    print(p)


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
