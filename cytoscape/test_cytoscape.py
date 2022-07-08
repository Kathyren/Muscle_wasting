import pytest

import Constants
import cytoscape as ct

cytoscape_small_data_nodes = Constants.cytoscape_small_data_nodes
cytoscape_small_data_edges = Constants.cytoscape_small_data_edges


def test_read_json():
    """
    This is a simple test to check if the lecture of the Cytoscape in in the right format
    :return:
    """
    js = ct.read_cytoscape_json(cytoscape_file='test.cyjs')
    assert js['elements'] == {'nodes': cytoscape_small_data_nodes,
                              'edges': cytoscape_small_data_edges}, f"The json file is in an unexpected format."


def test_format_cytoscape_json():
    cytoscape = Constants.cytoscape_small_string_json
    nodes, edges, relationsips = ct.format_cytoscape_json(cytoscape)
    assert len(edges) == len(
        relationsips), f"Given that the variable 'edges' is the metadata for the relationships, this must be the " \
                       f"exact same length "
    assert nodes == Constants.cytoscape_small_data_nodes
    assert edges == Constants.cytoscape_small_data_edges
    assert relationsips == Constants.cytoscape_small_relationships
