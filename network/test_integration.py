import Constants
import networkX as nx
import pytest
import main_pipeline as mp
import cytoscape as ct


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