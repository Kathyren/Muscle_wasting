import sys
sys.path.append('../')
sys.path.append('../network')
import pytest
from mirna_scoring.jupyter_functions import *
import mirna_scoring.mirna_impact as mi
import networkx as nx


def test_get_influences_df():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    df = mi.get_influences_df(graph)
    print (df)

    assert df

def test_get_minra_influence():


    pass

def test_get_impact_data():

    pass



def test_get_mirnas_similar_impact():
    pass

def test_cluster_mirnas():
    pass

def test_calculate_measurements():
    pass

def test_get_mirna_impact_pathway_rw():
    pass

