import sys
from operator import index

sys.path.append('../')
sys.path.append('../network')
import pytest
from mirna_scoring.jupyter_functions import *
import mirna_scoring.mirna_impact as mi
import networkx as nx
from mirna_scoring.walking_network import count_influence, evaluate_path_influence, get_influence


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
    """
    This function counts the keywords on the pathways
    :return:
    """
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    paths = {test_mir:[[test_mir, 'PKM'], [test_mir, 'ALDOA']]}
    path_pathway = mi.get_mirna_pathway_influence_df(network=graph, mirPaths=paths)
    print(path_pathway)
    assert path_pathway.shape == (1,10)
    pass

def test_get_mirna_pathway_influence_df_keywords():
    """
    This function counts the keywords on the pathways
    :return:
    """
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    paths = {test_mir:[[test_mir, 'PKM'], [test_mir, 'ALDOA']]}
    path_pathway = mi.get_mirna_pathway_influence_df(network=graph, mirPaths=paths, pathway_keywords=['ATP'])
    print(path_pathway)
    assert path_pathway.shape == (1,3)
    assert 'ATP' in path_pathway.columns
    pass


def test_get_mirna_influence():
    """
    This test is just to check that the influence goes on the paths for the only mir,
    :return:
    """
    test_mir = 'hsa-miR-122-5p'
    paths = {test_mir: [[test_mir, 'PKM'], [test_mir, 'ALDOA']]}
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    mirInfluence = mi.get_mirna_influence(mirPaths=paths, network=graph)
    print(mirInfluence)
    assert len(paths) == len(mirInfluence)
    assert len(mirInfluence[test_mir]) == 2
    m = {'hsa-miR-122-5p':
             [
                 {'dds_mf_m_file': 4.474424373018907,
                  'dds_mf_o_file': 2.9399779564138955,
                  'dds_mf_y_file': 3.150926539647954,
                  'dds_mo_F_file': 1.880962075599236,
                  'dds_mo_M_file': 2.3937433486756543,
                  'dds_mo_file': 2.077008111389381,
                  'dds_ym_F_file': 2.066059020071056,
                  'dds_ym_M_file': 0.7242056519795395,
                  'dds_ym_file': 1.5121426844660262,
                  'dds_yo_F_file': 3.3600807222800047,
                  'dds_yo_M_file': 1.5687836303279683,
                  'dds_yo_file': 4.794770137419901,
                  'pathways': ['GOBP_ATP_METABOLIC_PROCESS', 'GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY', 'GOBP_MONOCARBOXYLIC_ACID_METABOLIC_PROCESS', 'GOBP_NUCLEOBASE_CONTAINING_SMALL_MOLECULE_METABOLIC_PROCESS', 'GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS', 'GOBP_ORGANIC_ACID_METABOLIC_PROCESS', 'GOBP_ORGANONITROGEN_COMPOUND_BIOSYNTHETIC_PROCESS', 'GOBP_ORGANOPHOSPHATE_METABOLIC_PROCESS', 'GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS', 'GOBP_PYRUVATE_METABOLIC_PROCESS', 'GOBP_RIBOSE_PHOSPHATE_METABOLIC_PROCESS', 'GOBP_SMALL_MOLECULE_METABOLIC_PROCESS', 'GOCC_MITOCHONDRION', 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS', 'KEGG_PYRUVATE_METABOLISM'],
                  'pathways_svd': 456.7465652077005,
                  'weight': 67.97154212564476},
                 {'dds_mf_m_file': 4.08208246762406,
                  'dds_mf_o_file': 3.365215789627117,
                  'dds_mf_y_file': 4.025799822740248,
                  'dds_mo_F_file': 2.23564039293024,
                  'dds_mo_M_file': 2.1450823301001782,
                  'dds_mo_file': 2.356643793141844,
                  'dds_ym_F_file': 1.6945681463056437,
                  'dds_ym_M_file': 0.5774595982975826,
                  'dds_ym_file': 1.9840836691678088,
                  'dds_yo_F_file': 3.486283899282681,
                  'dds_yo_M_file': 1.484565112750286,
                  'dds_yo_file': 5.810481285595966,
                  'pathways': ['GOBP_ATP_BIOSYNTHETIC_PROCESS', 'GOBP_ATP_METABOLIC_PROCESS', 'GOBP_CARBOHYDRATE_DERIVATIVE_BIOSYNTHETIC_PROCESS', 'GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY', 'GOBP_MONOCARBOXYLIC_ACID_METABOLIC_PROCESS', 'GOBP_NUCLEOBASE_CONTAINING_SMALL_MOLECULE_METABOLIC_PROCESS', 'GOBP_NUCLEOSIDE_PHOSPHATE_BIOSYNTHETIC_PROCESS', 'GOBP_NUCLEOSIDE_TRIPHOSPHATE_BIOSYNTHETIC_PROCESS', 'GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS', 'GOBP_ORGANIC_ACID_METABOLIC_PROCESS', 'GOBP_ORGANONITROGEN_COMPOUND_BIOSYNTHETIC_PROCESS', 'GOBP_ORGANOPHOSPHATE_BIOSYNTHETIC_PROCESS', 'GOBP_ORGANOPHOSPHATE_METABOLIC_PROCESS', 'GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS', 'GOBP_PYRUVATE_METABOLIC_PROCESS', 'GOBP_RIBOSE_PHOSPHATE_BIOSYNTHETIC_PROCESS', 'GOBP_RIBOSE_PHOSPHATE_METABOLIC_PROCESS', 'GOBP_SMALL_MOLECULE_METABOLIC_PROCESS', 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS'],
                  'pathways_svd': 543.7033497448629,
                  'weight': 69.12395315378183}]}
    assert mirInfluence == m

def test_count_influence():
    test_mir = 'hsa-miR-122-5p'
    test_node = 'ALDOA'
    data_dict = {'weight': 0,
                 'pathways': [],
                 'pathways_svd': 0
                 }
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    node_data= graph.nodes[test_node]
    data_dict = count_influence(node_data=node_data, data_dict=data_dict, comparisons=None)
    print(data_dict)
    assert 'pathways' in data_dict
    assert 'pathways_svd' in data_dict
    assert 'weight' in data_dict
    pass
def test_count_influence_second_path():
    data_dict = {'dds_mf_m_file': 4.08208246762406, 'dds_mf_o_file': 3.365215789627117, 'dds_mf_y_file': 4.025799822740248,
                 'dds_mo_F_file': 2.23564039293024, 'dds_mo_M_file': 2.1450823301001782, 'dds_mo_file': 2.356643793141844,
                 'dds_ym_F_file': 1.6945681463056437, 'dds_ym_M_file': 0.5774595982975826, 'dds_ym_file': 1.9840836691678088,
                 'dds_yo_F_file': 3.486283899282681, 'dds_yo_M_file': 1.484565112750286, 'dds_yo_file': 5.810481285595966,
                 'pathways': ['GOBP_ATP_BIOSYNTHETIC_PROCESS', 'GOBP_ATP_METABOLIC_PROCESS'],
                 'pathways_svd': 543.7033497448629,
                 'weight': 18.123953153781827}
    test_node = 'ALDOA'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    node_data = graph.nodes[test_node]
    count_influence(node_data=node_data, data_dict=data_dict, comparisons=None)

    assert data_dict['pathways_svd'] == 2 * 543.7033497448629


def test_evaluate_path_influence():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    path = [test_mir, 'PKM', 'ALDOA', 'LDHA']
    eval = evaluate_path_influence(graph = graph, path=path)
    print(eval)
    assert eval['weight'] == 105.91091188330792

def test_get_influence():
    test_mir = 'hsa-miR-122-5p'
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    path = [test_mir, 'PKM', 'ALDOA', 'LDHA']
    paths = [[test_mir, 'PKM'], [test_mir, 'ALDOA']]
    influence = get_influence(graph=graph,paths=paths)
    print(influence)
    assert len(influence)== 2

    pass


def test_get_up_down_regulated():
    graph_pkl = '/home/karen/Documents/GitHub/Muscle_wasting/mirna_scoring/tests/test_files/ALDOA_LDHA_influence.pkl'
    graph = nx.read_gpickle(graph_pkl)
    u, d = mi.get_up_down_regulated(network=graph, condition='yo_file')
    print (u)
    assert u == {'ALDOA': 5.810481285595966, 'DNM1L': 2.8448790256653225,
                 'GAPDH': 6.441244380391132, 'LDHA': 6.986902151280292,
                 'PKM': 4.794770137419901}
    print(d)
    assert d == {'NT5E': -3.514083056569358}

def test_get_mirnas_similar_impact():
    pass
def test_cluster_mirnas():
    pass

def test_calculate_measurements():
    pass

def test_get_mirna_impact_pathway_rw():
    pass

def test_mirnas_influence_on_genes():
    network_path = "/home/karen/Documents/GitHub/Muscle_wasting/network/Networks_pkl/metadata_mirnas__Sarcopenia_relevant_normalize_cutoff_0.9.pkl"
    network = nx.read_gpickle(network_path)
    mirna = ['hsa-miR-16-5p']
    influence = mi.mirnas_influence_on_genes(network=network, miR_nodes=mirna)