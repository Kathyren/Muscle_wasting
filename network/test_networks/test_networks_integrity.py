import sys
sys.path.append('../')
#from .. import network_processing as nx
#from .. import main_pipeline as mp
import pytest
import cytoscape as ct
from network import main_pipeline as mp
from network import network_processing
import yaml
import pandas as pd

def test_test():
    assert True

def test_full_flow_genes_tf(monkeypatch):

    monkeypatch.setattr(
        mp, 'save_as_cjsn', lambda *args, **kwargs: None
     )
    name_n = "ALDOA_test"

    with open("network/settings/metadata.yml", 'r') as file:
        config_data = yaml.safe_load(file)

    file = config_data["oNetwork"]
    path_tissue_data = config_data['path_tissue_data']
    path_dds_data = config_data['path_DDS_data']
    pathway_file = config_data['path_pathway_file']
    tf_file = config_data['path_tf_file']


    import pandas as pd
    dds_df = pd.read_csv(path_dds_data, index_col=0).fillna(0)
    tissue_df = pd.read_csv(path_tissue_data, index_col=0).fillna(0)

    network, _ = mp.full_flow_genes_tf(cytoscape_network="network/test_networks/ALDOA_LDHA.cyjs",
                                    name= name_n,
                                    dds_df= dds_df,
                                    tissue_df=tissue_df,
                                    tf_file=tf_file,
                                    pathway_file=pathway_file,
                                    cutoff=0.5)
    print(network)
    assert len(network.nodes()) == 17, f"The test network ALDOA should have 6 genes to start and 11 mirnas."

def test_nodes_integrity(monkeypatch):
    graph = network_processing.load_graph("network/test_networks/original_ALDOA_test.pkl")
    monkeypatch.setattr(mp, "get_cytoscape_network",  lambda *args, **kwargs: graph )
    monkeypatch.setattr(
        mp, 'save_as_cjsn', lambda *args, **kwargs: None
     )
    name_n = "ALDOA_test"
    # monkeypatch for save_graph to do nothing
    monkeypatch.setattr(
        network_processing, 'save_graph', lambda *args, **kwargs: None
    )

    with open("network/settings/metadata.yml", 'r') as file:
        config_data = yaml.safe_load(file)

    file = config_data["oNetwork"]
    path_tissue_data = config_data['path_tissue_data']
    path_dds_data = config_data['path_DDS_data']
    pathway_file = config_data['path_pathway_file']
    tf_file = config_data['path_tf_file']


    import pandas as pd
    dds_df = pd.read_csv(path_dds_data, index_col=0).fillna(0)
    tissue_df = pd.read_csv(path_tissue_data, index_col=0).fillna(0)

    network, _ = mp.full_flow_genes_tf(cytoscape_network="network/test_networks/ALDOA_LDHA.cyjs",
                                    name= name_n,
                                    dds_df= dds_df,
                                    tissue_df=tissue_df,
                                    tf_file=tf_file,
                                    pathway_file=pathway_file,
                                    cutoff=0.5)
    print(network)
    assert len(network.nodes()) == 17, f"The test network ALDOA should have 6 genes to start and 11 mirnas."

def test_read_cytoscape_json(monkeypatch):
    """
    This test will cgheck that the nodes we start with are ok having the minimun properties.
    :param monkeypatch:
    :return:
    """
    cytoscape_json = ct.read_cytoscape_json(cytoscape_file="network/test_networks/ALDOA_LDHA.cyjs")
    nodes, edges, relationships = ct.format_cytoscape_json(cytoscape_json=cytoscape_json)

    print(nodes)
    # check that nodes is len 6
    assert len(nodes) == 6, f"The test network ALDOA should have 6 nodes."
    # check that edges is len 0
    assert len(edges) == 0, f"The test network ALDOA should have 0 edges."
    # check that relationships is len 0
    assert len(relationships) == 0, f"The test network ALDOA should have 0 relationships."
    # check that the nodes have the correct properties
    for node in nodes:
        assert 'data' in node, f"The node should have a data property"
        # Data must have name, id, shared_name, SUID and selected
        assert 'name' in node['data'], f"The node should have a name property"
        assert 'id' in node['data'], f"The node should have an id property"
        assert 'name' in node['data'], f"The node should have a name property"
        assert 'type' in node, f"The node should have a type property"
        assert node['type'] in ['gene', 'mirna'], f"The node should have a type property with value gene or mirna"

def test_get_cytoscape_network(monkeypatch):
    nodes = [{'data': {'SUID': 2000278, 'id': '2000278',
                       'id_original': '6071', 'name': 'NT5E',
                       'selected': False, 'shared_name': 'NT5E'},
              'position': {'x': 185.46056354998848, 'y': 296.23009851467975},
              'selected': False, 'source': 'Cytoscape', 'type': 'gene'},
             {'data': {'SUID': 2005417, 'id': '2005417',
                       'id_original': '2645', 'name': 'ALDOA',
                       'selected': False, 'shared_name': 'ALDOA'},
              'position': {'x': 57.46056354998848, 'y': 124.04747060021452},
              'selected': False, 'source': 'Cytoscape', 'type': 'gene'}]
    edges = []
    relationship= []

    # monkeypatch to read_cytoscape_json to retun nodes, edges and relationship
    monkeypatch.setattr("cytoscape.read_cytoscape_json", lambda *args, **kwargs: (nodes, edges, relationship))

    network = network_processing.create_graph_from_dictionaries(nodes=nodes, edges=edges, relationship=relationship)
    print (network)
    assert len(network.nodes()) == 2, f"The network should have 2 nodes"
    assert len(network.edges()) == 0, f"The network should have 0 edges"
    assert network.nodes['NT5E']['type'] == 'gene', f"The node NT5E should have type gene"
    assert network.nodes['ALDOA']['type'] == 'gene', f"The node ALDOA should have type gene"
