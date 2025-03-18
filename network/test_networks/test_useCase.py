import math
import network.network_processing as ntp
import numpy as np
def test_edge_weigth_values():
        network_path = "../../network/Networks_pkl/metadata_mirnas__Sarcopenia_relevant_normalize_cutoff_0.9.pkl"
        network = ntp.load_graph(network_path)
        ntp.weight_nodes(graph=network)
        ntp.weight_edges(graph=network, node_weight='weigh',
                         mir_enhancer=1)
        # for each edge in network get the weightScore, and get the min and max value
        min_value = 100000000
        max_value = 0
        nans = 0
        for edge in network.edges(data=True):
            values = edge[2]
            if 'weightScore' in values:
                if values['weightScore'] < min_value:
                    min_value = values['weightScore']
                if values['weightScore'] > max_value:
                    max_value = values['weightScore']
                if isinstance(values['weightScore'], float) and math.isnan(values['weightScore']):
                        nans += 1
            else:
                    print("No values")
        print(f"Min value: {min_value}")
        print(f"Max value: {max_value}")
        
        print(network)
def test_node_weight_values():
        network_path = "../../network/Networks_pkl/metadata_mirnas__Sarcopenia_relevant_normalize_cutoff_0.9.pkl"
        network = ntp.load_graph(network_path)
        ntp.weight_nodes(graph=network)
        min_value = 100000000
        max_value = 0
        nans = 0
        for node in network.nodes(data=True):
            values = node[1]['data']
            if 'weigh' in values:
                if values['weigh'] < min_value:
                    min_value = values['weigh']
                if values['weigh'] > max_value:
                    max_value = values['weigh']
                if isinstance(values['weigh'], float) and math.isnan(values['weigh']):
                        nans += 1
            else:
                    print("No values")
        print(f"Min value: {min_value}")
        print(f"Max value: {max_value}")
        print(f"Number of nans: {nans}")
        print(network)
        assert nans == 0, f"There should not be nans"
        assert min_value>=0, f"The min value should no be less than 0"

def test_mirna_attributes():
    network_path = "../../network/Networks_pkl/metadata_mirnas__Sarcopenia_relevant_normalize_cutoff_0.9.pkl"
    network = ntp.load_graph(network_path)
    # check all the nodes with type 'mirna' in the network
    mirnas = []
    for node in network.nodes(data=True):
        if 'type' in node[1]:            
                if node[1]['type'] == 'mirna':
                        mirnas.append(node)
        else:
            print("No type")
        # if node name is X, print 
        if node[0] == 'hsa-miR-16-5p':
            print(node)
    #mirnas = [node for node in network.nodes(data=True) if node[1]['type'] == 'mirna']
    print(mirnas)
    for mirna in mirnas:
        data = mirna[1]['data']
        print(data)
        assert data