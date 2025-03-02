import math

import pandas as pd
import os
from networkx import is_bipartite
from cytoscape.enum_network_sources import NetworkSource

import Constants
from network import network_processing
from network.main_pipeline import get_cytoscape_network
from network.network_processing import create_graph, draw_graph, load_graph, \
    create_graph_from_dictionaries, get_mirna_mrna_relationships, remove_nodes_low_centrality, convert_to_json, \
    save_graph, get_nodes_names, add_mirna_relationships, set_positions, get_mirna_tissue_edges, \
    add_tissue_relationship, add_organ_system_relationship, get_tissue_system_edges, extract_genes_from_pathways, \
    add_pathway_to_node, add_pathways_to_nodes, mark_TF_nodes_from_file, add_dds_to_node, add_DDS_data, random_walk

import pytest

genes = ['Cd320', 'Ndrg3', 'Aldoa', 'Bckdk', 'SLC7A1', 'ADAM17', 'NUMBL', 'FOXJ3', 'XPO6', 'AP3M2']
mirnas = ['mmu-miR-122-5p', 'mmu-miR-122-5p', 'mmu-miR-122-5p', 'mmu-miR-122-5p', 'hsa-miR-122-5p',
          'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p']
relationship = [('Cd320', 'mmu-miR-122-5p'), ('Ndrg3', 'mmu-miR-122-5p'), ('Aldoa', 'mmu-miR-122-5p'),
                ('Bckdk', 'mmu-miR-122-5p'), ('SLC7A1', 'hsa-miR-122-5p'), ('ADAM17', 'hsa-miR-122-5p'),
                ('NUMBL', 'hsa-miR-122-5p'), ('FOXJ3', 'hsa-miR-122-5p'), ('XPO6', 'hsa-miR-122-5p'),
                ('AP3M2', 'hsa-miR-122-5p')]
edges_metadata = []
cytoscape_nodes = [{
    "data": {
        "id": "166",
        "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
        "IntAct_AC": "EBI-10206780",
        "IntAct_Mutation": False,
        "IntAct_Taxon_Id": "9606",
        "IntAct_Preferred_Id": "P35523",
        "IntAct_Preferred_Id_Database": "uniprotkb",
        "IntAct_Species": "Homo sapiens",
        "IntAct_Type_MI_identifier": "MI:0326",
        "shared_name": "CLCN1",
        "IntAct_Type": "protein",
        "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
        "IntAct_Identifier_Identifiers": ["EBI-10206783", "EBI-10206784", "EBI-30414706", "EBI-10206780",
                                          "EBI-10206782"],
        "IntAct_Description": "Chloride channel protein 1",
        "name": "CLCN1",
        "SUID": 166,
        "selected": False,
        "IntAct_Feature_Features": ["EBI-25146010", "EBI-25065318", "EBI-24991613", "EBI-23932890", "EBI-23509792",
                                    "EBI-23122996"]
    },
    "position": {
        "x": -48.68993826842569,
        "y": 35.43257138655283
    },
    "selected": False}, {
    "data": {
        "id": "164",
        "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
        "IntAct_AC": "EBI-7062247",
        "IntAct_Mutation": False,
        "IntAct_Taxon_Id": "9606",
        "IntAct_Preferred_Id": "Q9UHD4",
        "IntAct_Preferred_Id_Database": "uniprotkb",
        "IntAct_Species": "Homo sapiens",
        "IntAct_Type_MI_identifier": "MI:0326",
        "shared_name": "CIDEB",
        "IntAct_Type": "protein",
        "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
        "IntAct_Identifier_Identifiers": ["EBI-30215058", "EBI-30215059", "EBI-30215060", "EBI-30215061",
                                          "EBI-30215062", "EBI-30215063", "EBI-7062251", "EBI-7062252",
                                          "EBI-7062253", "EBI-7062247", "EBI-7062250"],
        "IntAct_Description": "Cell death activator CIDE-B",
        "name": "CIDEB",
        "SUID": 164,
        "selected": True,
        "IntAct_Feature_Features": ["EBI-25146017", "EBI-23932897"]
    },
    "position": {
        "x": -54.39876059290216,
        "y": 145.18372932861945
    },
    "selected": True}, {
    "data": {
        "id": "162",
        "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
        "IntAct_AC": "EBI-12859340",
        "IntAct_Mutation": False,
        "IntAct_Taxon_Id": "9606",
        "IntAct_Preferred_Id": "Q9NQX1-2",
        "IntAct_Preferred_Id_Database": "uniprotkb",
        "IntAct_Species": "Homo sapiens",
        "IntAct_Type_MI_identifier": "MI:0326",
        "shared_name": "PRDM5",
        "IntAct_Type": "protein",
        "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
        "IntAct_Identifier_Identifiers": ["EBI-30096885", "EBI-12859340", "EBI-12859341"],
        "IntAct_Description": "PR domain zinc finger protein 5",
        "name": "PRDM5",
        "SUID": 162,
        "selected": False,
        "IntAct_Feature_Features": ["EBI-25065325", "EBI-23509799"]
    },
    "position": {
        "x": -201.89180646623885,
        "y": 33.63696071079798
    },
    "selected": False
}]

graph = load_graph("network/small_graph.pkl")



def test_create_network():
    graph = create_graph(genes=genes, mirnas=mirnas, relationsip=relationship)
    assert len(graph.nodes()) == 12
    assert len(graph.edges()) == 10
    draw_graph(graph)


def test_draw():
    graph = create_graph(genes=genes, mirnas=mirnas, relationsip=relationship)
    draw_graph(graph)
    pass


def test_is_bipartite():
    graph = load_graph("small_graph.pkl")
    x = is_bipartite(graph)


def test_remove_nodes_low_centrality():
    graph = load_graph("big_graph.pkl")
    n_nodes = graph.number_of_nodes()
    graph_p = remove_nodes_low_centrality(graph)
    n_nodes_2 = graph_p.number_of_nodes()
    assert n_nodes >= n_nodes_2


def test_create_graph_from_dictioanries():
    nodes = Constants.cytoscape_small_data_nodes
    edges = Constants.cytoscape_small_data_edges
    ppi = Constants.cytoscape_small_string_json
    graph = create_graph_from_dictionaries(nodes=nodes, edges=edges, relationship=relationship)
    print(list(graph.nodes))
    print(list(graph.edges))


def test_get_mirna_mrna_relationships():
    relationship, mirnas, scores = get_mirna_mrna_relationships(genes=genes)
    assert type(relationship) is list
    assert type(relationship[0]) is tuple
    assert relationship[0][0] in genes
    assert relationship[0][1] in mirnas
    assert type(scores) is list
    assert scores[0] > 0


def test_convert_to_json():
    graph = load_graph("big_graph.pkl")
    dict_graph = convert_to_json(graph)

    assert not (dict_graph.get(
        'elements') is None), f"Json graph do not have essential key \"element\", the keys are: {dict_graph.keys} "
    assert not (dict_graph['elements'].get('nodes') is None), f"Json graph do not have essential key \"nodes\", " \
                                                              f"the keys are: {dict_graph['elements'].keys} "
    assert not (dict_graph['elements'].get('edges') is None), f"Json graph do not have essential key \"edges\", " \
                                                              f"the keys are: {dict_graph['elements'].keys} "


def test_convert_network_to_json():
    nodes = Constants.cytoscape_small_data_nodes
    edges = Constants.cytoscape_small_data_edges
    ppi = Constants.cytoscape_small_string_json
    graph = create_graph_from_dictionaries(nodes=nodes, edges=edges, relationship=relationship)
    dict_graph = convert_to_json(graph)

    assert not (dict_graph.get(
        'elements') is None), f"Json graph do not have essential key \"element\", the keys are: {dict_graph.keys} "
    assert not (dict_graph['elements'].get('nodes') is None), f"Json graph do not have essential key \"nodes\", " \
                                                              f"the keys are: {dict_graph['elements'].keys} "
    assert not (dict_graph['elements'].get('edges') is None), f"Json graph do not have essential key \"edges\", " \
                                                              f"the keys are: {dict_graph['elements'].keys} "


def test_get_nodes_names():
    nodes = Constants.cytoscape_small_data_nodes
    node_keys = Constants.cytoscape_small_network_node_names
    edges = Constants.cytoscape_small_data_edges
    relationship = Constants.cytoscape_small_relationships
    graph = create_graph_from_dictionaries(nodes=nodes, edges=edges, relationship=relationship)
    node_names = get_nodes_names(graph)
    assert set(node_keys) == set(node_names)


def test_add_mirnas(monkeypatch):
    def fake_relationships(*args, **kwargs):
        genes = Constants.cytoscape_small_network_node_names
        fake_relationship = [(genes[0], 'hsa-miR-111-5p'), (genes[1], 'hsa-miR-111-5p'),
                             (genes[1], 'hsa-miR-122-3p'), (genes[1], 'hsa-miR-130-5p'), (genes[1], 'hsa-miR-122-5p'),
                             (genes[0], 'hsa-miR-122-5p')]
        return fake_relationship, ['hsa-miR-111-5p', 'hsa-miR-122-3p',
                                   'hsa-miR-130-5p', 'hsa-miR-122-5p'], [1] * len(relationship)

    monkeypatch.setattr(
        networkX, "get_mirna_mrna_relationships", fake_relationships
    )
    node_count_1 = graph.number_of_nodes()
    add_mirna_relationships(graph)
    node_count_2 = graph.number_of_nodes()
    assert node_count_2 > node_count_1, f'No added nodes?'
    # save_graph(graph,'test_w_mirnas.pkl')


def test_set_positions():
    x = set_positions(graph)
    a = Constants.cytoscape_small_network_node_names[0]
    b = Constants.cytoscape_small_network_node_names[1]
    pos_a = [x[a][0], x[a][1]]
    pos_b = [x[b][0], x[b][1]]
    assert graph.nodes[a]['position']['x'] == pos_a[0], f"The value for {a} in x is {pos_a[0]}," \
                                                        f" {graph.nodes[a]['position']['x']} found"
    assert graph.nodes[a]['position']['y'] == pos_a[1], f"The value for {a} in y is {pos_a[1]}," \
                                                        f" {graph.nodes[a]['position']['y']} found"
    assert graph.nodes[b]['position']['x'] == pos_b[0], f"The value for {b} in x is {pos_a[0]}," \
                                                        f" {graph.nodes[b]['position']['x']} found"
    assert graph.nodes[b]['position']['y'] == pos_b[1], f"The value for {b} in y is {pos_a[1]}," \
                                                        f" {graph.nodes[b]['position']['y']} found"


def test_get_mirna_tissue_edges_random_mirnas():
    random_mirnas = ['hsa-miR-181a-2-3p', 'hsa-miR-519b-5p', 'hsa-miR-629-5p', 'hsa-miR-6721-5p', 'hsa-miR-519a-3p']
    edges, organs, source = get_mirna_tissue_edges(random_mirnas)
    assert len(edges) > 1
    assert source[0] > 0


def test_add_tissue_relationship(monkeypatch):
    def fake_relationships(*args, **kwargs):
        fake_relationship = [('muscle', 'hsa-miR-111-5p'), ('brain', 'hsa-miR-111-5p'),
                             ('muscle', 'hsa-miR-122-3p'), ('vein', 'hsa-miR-130-5p'), ('bone', 'hsa-miR-122-5p'),
                             ('vein', 'hsa-miR-122-5p')]
        return fake_relationship, ['muscle', 'vein',
                                   'bone', 'brain'], [1] * len(relationship)

    monkeypatch.setattr(
        networkX, "get_mirna_tissue_edges", fake_relationships
    )
    graph = load_graph("test_w_mirnas.pkl")
    node_count_1 = graph.number_of_nodes()
    add_tissue_relationship(graph)
    node_count_2 = graph.number_of_nodes()
    assert node_count_2 >= node_count_1 + 4, f'No added nodes?'
    # (graph,'test_w_mirnas_organs.pkl')
    cc_it = networkx.all_neighbors(graph=graph, node='muscle')
    cc = len(list(cc_it))
    # save_graph(graph, "test_w_mirnas_organs.pkl")
    assert cc == 2, f"The fake relationships is designed to give 2 neighbors for muscle." \
                    f" Maybe mocked graph already had muscles? "


def test_add_organ_system_relationship(monkeypatch):
    def fake_relationships(*args, **kwargs):
        fake_relationship = [('muscle', 'lymphatic'), ('brain', 'lymphatic'),
                             ('musculoskeletal', 'muscle'), ('vein', 'lymphatic'), ('bone', 'musculoskeletal'),
                             ('vein', 'cardiovascular')]
        return fake_relationship, ['lymphatic', 'musculoskeletal',
                                   'cardiovascular']

    monkeypatch.setattr(
        networkX, "get_tissue_system_edges", fake_relationships
    )
    graph = load_graph("test_w_mirnas_organs.pkl")
    node_count_1 = graph.number_of_nodes()
    add_organ_system_relationship(graph)
    node_count_2 = graph.number_of_nodes()
    assert node_count_2 > node_count_1, f'No added nodes?'

    # save_graph(graph, "test_w_mirnas_organs_n_systems.pkl")
    cc_it = networkx.all_neighbors(graph=graph, node='lymphatic')
    cc = len(list(cc_it))
    assert cc == 3, f"The fake relationships is designed to give 2 neighbors for muscle." \
                    f" Maybe mocked graph already had muscles? "


def test_get_tissue_system_edges():
    organs = ['muscle', 'brain']
    edges, systems = get_tissue_system_edges(organs=organs)
    assert len(edges) > 1
    assert edges[0][0] in organs


def test_get_network():
    network = "graph1_Selected_genes2.cyjs"
    the_network = get_cytoscape_network(network)
    nodes = the_network.nodes


def test_extract_nodes_from_pathways():
    result = extract_genes_from_pathways('test_enr_pathways.csv')
    print(result)
    assert len(result) == 98, "not all the genes are registered here"
    assert 'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA' in result['ATP5F1A']
    assert 'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA' in result['COX17']


def test_add_pathway_to_node(monkeypatch):
    graph = load_graph("network/graph0_tf_network_cutoff_0.pkl")
    test_node = dict(graph.nodes(data=True))

    node_name = 'BGLAP'
    test_node = test_node[node_name]
    pathways = ['P1', 'P2']

    add_pathway_to_node(graph, node_name, pathways)

    assert test_node['data']['pathways'] == pathways


def test_add_dds_to_node():
    graph = load_graph("network/graph0_tf_network_cutoff_0.pkl")
    row_dict = {'yo': -0.266, 'ym': None, 'mo': None, 'ml_c': None, 'ml_s': None}
    node = 'CFH'
    add_dds_to_node(graph=graph, node_name=node, dds=dict(row_dict))
    test_nodes = dict(graph.nodes(data=True))
    dds = test_nodes[node]
    assert dds['data']['yo'] == -0.266
    # check if the property metadata is in the node
    assert 'data' in dds
    assert 'metadata' in dds['data']
    assert 'dds' in dds['data']['metadata']
    assert dds['data']['metadata']['dds'] == row_dict
    pass


def test_add_dds_to_node_non_existent():
    graph = load_graph("graph0_tf_network_cutoff_0.pkl")
    row_dict = {'yo': -0.266, 'ym': None, 'mo': None, 'ml_c': None, 'ml_s': None}
    node = 'HHH'
    add_dds_to_node(graph=graph, node_name=node, dds=dict(row_dict))
    test_nodes = dict(graph.nodes(data=True))
    assert node not in test_nodes


def test_extract_genes_from_pathways():
    pathway_file = 'network/test_networks/dummy_data/pathway_file.csv'
    pathway_df = pd.read_csv(pathway_file)
    feature_dict = network_processing.extract_genes_from_pathways(pathway_df=pathway_df , id_feature = 'genes')
    assert len(feature_dict['MYC']) == 2
    assert len(feature_dict['NT5E']) == 1
    assert feature_dict['NT5E'][0] == 'GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE'
    pass

def test_add_pathways_to_nodes():
    graph = load_graph("network/graph0_tf_network_cutoff_0.pkl")
    pathway_file = 'network/test_networks/dummy_data/pathway_file.csv'
    add_pathways_to_nodes(graph, pathway_file)
    test_nodes = dict(graph.nodes(data=True))
    #t2 = dict(g2.nodes(data=True))
    print(test_nodes)
    assert len(test_nodes['MYC']['data']['pathways']) == 2
    # assert metadata in data
    assert 'metadata' in test_nodes['MYC']['data']
    assert 'pathways' in test_nodes['MYC']['data']['metadata']
    assert len(test_nodes['MYC']['data']['metadata']['pathways']) == 2
    

def test_get_SVD_to_node():
    pathway_df = pd.read_csv('network/test_networks/dummy_data/pathway_file.csv')
    pathway_list = ['GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE', 
                    'GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE']
    pathway_df = pathway_df.drop(columns=['genes'])
    pv, svd, cv = network_processing.get_SVD_pathways(pathway_df=pathway_df, pathways_list=pathway_list)
    print(pv)
    print(svd)
    assert pv == [-1]
    assert math.ceil(svd)==394


def test_add_tf_to_node():
    graph = load_graph("network/graph0_tf_network_cutoff_0.pkl")
    row_dict = {'yo': -0.266, 'ym': None, 'mo': None, 'ml_c': None, 'ml_s': None}
    node = 'CFH'
    network_processing.add_tf_to_node(graph=graph, node_name=node, tf=dict(row_dict))
    test_nodes = dict(graph.nodes(data=True))
    tf = test_nodes[node]
    # check if the property metadata is in the node
    assert 'data' in tf
    assert 'metadata' in tf['data']
    assert 'tf' in tf['data']['metadata']
    assert tf['data']['metadata']['tf'] == row_dict


def test_add_DDS_data():
    graph = load_graph("graph0_tf_network_cutoff_0.pkl")

    ddf_dict_list = [{'gene': 'CFH', 'yo': -0.2, 'ym': 3, 'mo': None, 'ml_c': None, 'ml_s': -0.5},
                     {'gene': 'SEMA3F', 'yo': -4.10, 'ym': None, 'mo': -3, 'ml_c': None, 'ml_s': None},
                     {'gene': 'CFTR', 'yo': -3.34, 'ym': None, 'mo': -4, 'ml_c': None, 'ml_s': 0.5}]
    df = pd.DataFrame.from_dict(ddf_dict_list)
    df = df.set_index('gene')
    add_DDS_data(graph=graph, dds_df=df)
    test_nodes = dict(graph.nodes(data=True))
    print(test_nodes)
    assert test_nodes['CFH']['data']['yo'] == -0.2
    assert test_nodes['CFTR']['data']['ml_s'] == 0.5
    pass


def test_mark_TF_nodes_from_file():
    graph = load_graph("graph0_tf_network_cutoff_0.pkl")
    tf_file = 'test_tf_act.csv'
    mark_TF_nodes_from_file(graph, tf_file)

    nodes = dict(graph.nodes(data=True))
    test_node = nodes['RELA']
    assert test_node['data']['node_type'] == 'TF'


def test_random_walk_no_neig():
    """
    Test random walk of a node with no out nodes
    :return:
    """
    graph = load_graph("graph0_tf_network_cutoff_0.pkl")
    gene = 'CFTR'

    paths = random_walk(graph=graph, node_name=gene, distance=5, sample_size=1)
    print(paths)
    assert len(paths) == 1
    a = paths[0]
    assert a == []

def test_random_walk():
    """
    Test random walk of a node
    :return:
    """
    graph = load_graph("mirnas_tf_network_cutoff_0.pkl")
    gene = 'hsa-miR-21-5p'
    sample_size = 100
    n_distance = 10
    paths = random_walk(graph=graph, node_name=gene, distance=n_distance, sample_size=sample_size)
    print(paths)
    assert len(paths) == sample_size
    a = paths[0]
    print(a)
    x = list(filter(lambda i: len(i) > 3, paths))
    print(x)
    assert len(a)>=1
    assert len(a)<=n_distance

def test_random_walk_leave_pathway():
    """
    Test random walk of a node
    :return:
    """
    graph = load_graph("mirnas_tf_network_cutoff_0.pkl")
    gene = 'hsa-miR-584-5p'
    sample_size = 100
    n_distance = 10
    paths = random_walk(graph=graph, node_name=gene, distance=n_distance, sample_size=sample_size)
    print(paths)
    assert len(paths) == sample_size
    a = paths[0]
    print(a)
    x = list(filter(lambda i: len(i) > 3, paths))
    print(x)
    assert len(a)>=1
    assert len(a)<=n_distance

