import pandas as pd
import pytest

import mirna_scoring.mirna_impact as mi
import networkx as nx
graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
graph = nx.read_gpickle(graph_pkl)
my_network = mi.mirna_network(network=graph)

def test_initialize_object_miR_nodes():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mi_network = mi.mirna_network(network=graph)
    miR_nodes = mi_network.miR_nodes
    assert miR_nodes == [test_mir]

def test_initialize_object_impact_df():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mi_network = mi.mirna_network(network=graph)
    influence_df = mi_network.influence_df
    print(influence_df)
    nodes = graph.nodes
    ALDOA = nodes['ALDOA']['data']['influence'][test_mir]
    PKM = nodes['PKM']['data']['influence'][test_mir]
    assert ALDOA == [-1]
    assert PKM == [-1]
    assert 'influence' not in nodes['NT5E']['data'], f"No mirna influences NT5E"
def test_initialize_object_influence_sum_df():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mi_network = mi.mirna_network(network=graph)
    influence_sum_df = mi_network.influence_sum_df
    print(influence_sum_df)
    nodes = graph.nodes
    influence = influence_sum_df['ALDOA']['hsa-miR-122-5p']
    assert influence == -1
    influence = influence_sum_df['PKM']['hsa-miR-122-5p']
    assert influence == -1

def test_up_down_regulated():
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_end.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mi_network = mi.mirna_network(network=graph)
    condition = 'yo_file'
    mi_network.set_up_down_regulated(condition=condition)
    up_regulated = mi_network.upRegulated
    print(up_regulated)
    assert up_regulated == {'yo_file': {'ALDOA': 5.810481285595966, 'DNM1L': 2.8448790256653225,
                                        'GAPDH': 6.441244380391132, 'LDHA': 6.986902151280292,
                                        'PKM': 4.794770137419901}}
    down_regulated = mi_network.downRegulated
    assert down_regulated == {'yo_file': {'NT5E': -3.514083056569358}}

    assert up_regulated[condition] == mi_network.get_up_regulated(condition=condition)
    assert down_regulated[condition] == mi_network.get_down_regulated(condition=condition)

def test_set_all_conditions_up_down_regulated():
    my_network.set_all_conditions_up_down_regulated()
    up_reg = my_network.get_up_regulated('yo_file')
    assert up_reg == {'ALDOA': 5.810481285595966, 'DNM1L': 2.8448790256653225,
                            'GAPDH': 6.441244380391132, 'LDHA': 6.986902151280292,
                            'PKM': 4.794770137419901}

def test_get_random_walk_pathway_influence():
    mir_pathway_influence_df = my_network.get_random_walk_pathway_influence()
    print(mir_pathway_influence_df)
    test_mir = 'hsa-miR-122-5p'
    miR_122 = mir_pathway_influence_df.loc[test_mir]
    assert miR_122['ATP']==3
    assert miR_122['RESPIRAT']==0

def test_calculate_measurements():
    measurements = my_network.calculate_measurements(comparisons=['yo_file', 'mo_file', 'ym_file'])
    print(measurements)
    for measure, measurement in measurements.items():
        print(measurement.shape)
        assert measurement.shape==(1,2), f"There is one microRNA and 2 genes, all dataframes should be the same outlet"


def test_all_scorings():
    my_network.set_all_conditions_up_down_regulated()
    up = my_network.calculate_up_regulated_dds_score()
    my_network.set_up_regulated_dds_score()
    down = my_network.calculate_down_regulated_dds_score()
    my_network.set_down_regulated_dds_score()
    pathways_sdv = my_network.calculate_pathway_score()
    my_network.set_pathway_score()
    pathway_count = my_network.get_random_walk_pathway_influence(pathway_keywords=['ATP'])
    all = my_network.get_scores()
    assert pd.DataFrame(up).equals(all['activators'])
    assert pd.DataFrame(down).equals(all['inhibitors'])
    print (up)
