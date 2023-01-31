import Constants
import networkX as nx
import pytest
import main_pipeline as mp
import cytoscape as ct
from network import main_pipeline


def test_get_readable_output(monkeypatch):
    """
    01
    This is actually a manual test to see if whatever we are getting is readable by cytoscape
    :return:
    """

    def fake_relationships(*args, **kwargs):
        genes = Constants.cytoscape_small_network_node_names
        fake_relationship = [(genes[0], 'hsa-miR-111-5p'), (genes[1], 'hsa-miR-111-5p'),
                             (genes[1], 'hsa-miR-122-3p'), (genes[1], 'hsa-miR-130-5p'), (genes[1], 'hsa-miR-122-5p'),
                             (genes[0], 'hsa-miR-122-5p')]
        return fake_relationship, ['hsa-miR-111-5p', 'hsa-miR-122-3p',
                                   'hsa-miR-130-5p', 'hsa-miR-122-5p'], [1] * len(fake_relationship)

    monkeypatch.setattr(
        nx, "get_mirna_mrna_relationships", fake_relationships
    )
    graph = nx.load_graph("small_graph.pkl")
    nx.add_mirna_relationships(graph)
    nx.set_positions(graph)
    mp.save_as_cjsn(graph, f'test_integration_01.cyjs')
    graph2 = ct.read_cytoscape_json("test_integration_01.cyjs")
    assert graph2


def test_main(monkeypatch):
    def fake_relationships(*args, **kwargs):
        genes = Constants.cytoscape_small_network_node_names
        fake_relationship = [(genes[0], 'hsa-miR-111-5p'), (genes[1], 'hsa-miR-111-5p'),
                             (genes[1], 'hsa-miR-122-3p'), (genes[1], 'hsa-miR-130-5p'), (genes[1], 'hsa-miR-122-5p'),
                             (genes[0], 'hsa-miR-122-5p')]
        return fake_relationship, ['hsa-miR-111-5p', 'hsa-miR-122-3p',
                                   'hsa-miR-130-5p', 'hsa-miR-122-5p'], [1] * len(fake_relationship)

    monkeypatch.setattr(
        nx, "get_mirna_mrna_relationships", fake_relationships
    )
    mp.add_mirnas_n_select("test_integration_01.cyjs", "dryrun_small_w_mirnas", use_prefix=False)



def test_add_mirnas_n_tissues():
    cytoscape_network = "STRING_Magagnes2009.cyjs"

    mp.add_mirnas_n_tissues(cytoscape_network=cytoscape_network, name="test_mirnas_tissues_evaluating")



def test_remove():
    cutoff = 0.85
    the_network = nx.load_graph(f"graph1_Selected_genes.pkl")
    save_as_cjsn(the_network, f'check_Selected_genes.cyjs')
    network = nx.remove_nodes_low_centrality(graph=the_network, cutoff=cutoff)
    nx.save_graph(network, f"graph3_Selected_genes.pkl")
    save_as_cjsn(network, f'graph1_Selected_genes.cyjs')
    network2= get_cytoscape_network("graph1_Selected_genes.cyjs")
    pass

def test_check_sanity():
    the_network = nx.load_graph(f"graph1_Selected_genes.pkl")
    nodes = the_network._node
    edges = the_network.edges()

    for edge in edges:
        print(edge)
        if edge[0] in nodes and edge[1] in nodes:
            pass
        else:
            print(f"edge {edge[0]}-{edge[1]} ")
    #nx.save_graph(the_network, f"test.pkl")
    mp.save_as_cjsn(the_network, f'test.cyjs')
    network2= mp.get_cytoscape_network("test.cyjs")
    x = network2.nodes

    pass

def test_node_integration(monkeypatch):
    """
    This test is to check if the nodes are keeping the same values after being treated
    :return:
    """

    def fake_relationships(*args, **kwargs):

        genes = ['EPB41L3', 'BDKRB2', 'ATF3', 'Ywhah']
        fake_relationship = [(genes[0], 'hsa-miR-111-5p'), (genes[1], 'hsa-miR-111-5p'),
                             (genes[1], 'hsa-miR-122-3p'), (genes[1], 'hsa-miR-130-5p'), (genes[1], 'hsa-miR-122-5p'),
                             (genes[0], 'hsa-miR-122-5p')]
        return fake_relationship, ['hsa-miR-111-5p', 'hsa-miR-122-3p',
                                   'hsa-miR-130-5p', 'hsa-miR-122-5p'], [1] * len(fake_relationship)

    #monkeypatch.setattr(
    #    nx, "get_mirna_mrna_relationships", fake_relationships
    #)
    #monkeypatch.setattr(
    #    nx, 'save_graph', lambda x, y: None
    #)
    monkeypatch.setattr(
        nx, 'remove_nodes_low_centrality', lambda graph, cutoff: graph
    )
    #monkeypatch.setattr(
    #    main_pipeline, 'save_as_cjsn', lambda *args, **kwargs: None
    #)
    monkeypatch.setattr(
        nx, 'set_positions', lambda *args, **kwargs: None
    )

    file_n = "/home/karen/Documents/GitHub/Muscle_wasting/cytoscape/Diff_express_genes.cyjs"
    name_n = "Selected_genes_testing"
    network = mp.add_mirnas_n_select(file_n, name_n, cutoff=0.85)
    nodes = network.nodes()
    edges = network.edges()
    network2 = mp.get_cytoscape_network(f'graph1_{name_n}.cyjs')
