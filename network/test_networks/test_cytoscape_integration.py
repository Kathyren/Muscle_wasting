import sys
sys.path.append('../')
import pytest
import cytoscape as ct
from network import main_pipeline as mp
from network import network_processing

import yaml
import pandas as pd



def test_format_cytoscape_json(monkeypatch):
    """
    This test will cgheck that the nodes we start with are ok having the minimun properties.
    :param monkeypatch:
    :return:
    """
    network = mp.get_cytoscape_network(file_name="network/test_networks/ALDOA_LDHA.cyjs")
    # Check the number of nodes be 6
    assert len(network.nodes()) == 6, f"The test network ALDOA should have 6 genes to start."
    # check that the nodes have the attributes data, position and selected
    for node in network.nodes():
        assert "data" in network.nodes[node], f"The node {node} should have the attribute data"
        assert "position" in network.nodes[node], f"The node {node} should have the attribute position"
        assert "selected" in network.nodes[node], f"The node {node} should have the attribute selected"
        assert "type" in network.nodes[node], f"The node {node} should have the attribute type"


