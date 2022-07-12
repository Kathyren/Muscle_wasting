from networkx import is_bipartite

import Constants
from network import networkX
from network.main_pipeline import get_cytoscape_network
from network.networkX import create_graph, draw_graph, load_graph, \
    create_graph_from_dictionaries, get_mirna_mrna_relationships, remove_nodes_low_centrality, convert_to_json, \
    save_graph, get_nodes_names, add_mirna_relationships
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
graph = load_graph("small_graph.pkl")


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
    relationship = get_mirna_mrna_relationships(genes=genes)
    assert type(relationship) is list
    assert type(relationship[0]) is tuple
    assert relationship[0][0] in genes


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
        return fake_relationship, ['hsa-miR-111-5p', 'hsa-miR-122-3p', 'hsa-miR-130-5p', 'hsa-miR-122-5p']

    monkeypatch.setattr(
        networkX, "get_mirna_mrna_relationships", fake_relationships
    )
    node_count_1 = graph.number_of_nodes()
    add_mirna_relationships(graph)
    node_count_2 = graph.number_of_nodes()
    assert node_count_2 > node_count_1, f'No added nodes?'
