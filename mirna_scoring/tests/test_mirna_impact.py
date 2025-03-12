import sys
from operator import index

sys.path.append('../')
sys.path.append('../network')
import pytest
from mirna_scoring.jupyter_functions import *
import mirna_scoring.mirna_impact as mi
import networkx as nx


def test_get_influences_df():
    """
    This test is to check in a simple case, the influcne DF is correct. In ALDOA_LDHA_influence the nodes
    already have the influence attribute. Only ALDOA and PKM have edges.
    :return:
    """
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    df = mi.get_influences_df(graph)
    assert len(df) == 1, "There is only one mirna affecting"
    assert df['ALDOA'][test_mir]==[-1]

def test_get_minra_influence():
    influence_df = pd.DataFrame(columns=['ALDOA', 'PKM'], data=[[[-1], [-1]]], index=['hsa-miR-122-5p'])
    test_mir = 'hsa-miR-122-5p'
    mir_influence = get_minra_influence(test_mir=test_mir, df = influence_df)
    print(mir_influence)
    tr = influence_df.T
    assert len(mir_influence) == 2
    assert mir_influence.equals(tr)


def test_get_impact_data_trivial_case():
    influence_df = pd.DataFrame(columns=['ALDOA', 'PKM'], data=[[[-1], [-1]]], index=['hsa-miR-122-5p'])
    int_df =  get_impact_data(influence_df)
    print(int_df)
    assert int_df.shape == influence_df.shape

def test_get_impact_data():
    influence_df = pd.DataFrame(columns=['ALDOA', 'PKM'], data=[[[-1, -1, 1], [-1, 1]]], index=['hsa-miR-122-5p'])
    int_df =  get_impact_data(influence_df)
    print(int_df)
    assert int_df.shape == influence_df.shape
    influence = int_df['ALDOA']['hsa-miR-122-5p']
    assert influence == -1
    influence = int_df['PKM']['hsa-miR-122-5p']
    assert influence == 0

def test_get_paths():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    paths = mi.get_paths(network=graph, nodes_start=[test_mir], steps=5, sample_size=40)
    print(paths)
    mir_paths = paths[test_mir]
    assert len(mir_paths) == 2
    assert [test_mir, 'PKM'] in mir_paths
    assert [test_mir, 'ALDOA'] in mir_paths

    pass



def test_get_mirna_pathway_influence_df():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    paths = {test_mir:[[test_mir, 'PKM'], [test_mir, 'ALDOA']]}
    path_pathway = mi.get_mirna_pathway_influence_df(network=graph, mirPaths=paths)
    print(path_pathway)
    assert path_pathway.shape == (1,10)
    pass


def test_get_mirna_influence():
    test_mir = 'hsa-miR-122-5p'
    paths = {test_mir: [[test_mir, 'PKM'], [test_mir, 'ALDOA']]}
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mirInfluence = mi.get_mirna_influence(mirPaths=paths, network=graph)
    print(mirInfluence)
    pass




def test_get_mirnas_similar_impact():
    pass
def test_cluster_mirnas():
    pass

def test_calculate_measurements():
    pass

def test_get_mirna_impact_pathway_rw():
    pass

