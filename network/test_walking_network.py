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
