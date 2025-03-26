import pytest
import mirkitten.create_network as cn
import pandas as pd
import os

def test_initialize_network():
    dds_files_path = 'dummy_data/dummy_dds.yml'
    # print pwd
    pwd = os.getcwd()
    print(pwd)
    network = cn.new_network(dds_files_path, save_name='dummy_results/network_test.cyjs')
    assert network.only_DE == False
    assert network.de_threshold == None
    assert network.pvalue == 0.05
    assert network.interest == 'stat'
    assert network.species == 'human'
    assert network.source_target_df == None
def test_get_nodes():
    dds_files_path = 'dummy_data/dummy_dds.yml'
    network = cn.new_network(dds_files_path, save_name='dummy_results/network_test.cyjs')
    network.collect_genes()
    nodes= network.nodes
    print(nodes)
    assert 'TSPAN6' in nodes
    assert 'TNMD' in nodes
    assert 'DPM1' in nodes
    assert 'SCYL3' in nodes
def test_get_network():
    dds_files_path = 'dummy_data/dummy_dds.yml'
    network = cn.new_network(dds_files_path, save_name='dummy_results/network_test.cyjs')
    network.collect_genes()
    network.get_network()
    source_target_df = network.source_target_df
    print(source_target_df)
    assert 'source' in source_target_df.columns
    assert 'target' in source_target_df.columns
    assert 'weight' in source_target_df.columns
    assert 'PMID' in source_target_df.columns
    assert len(source_target_df)==5
def test_cytoscape_network():
    dds_files_path = 'dummy_data/dummy_dds.yml'
    network = cn.new_network(dds_files_path,
                             save_name='dummy_results/network_test.cyjs')
    network.collect_genes()
    network.get_network()
    cytoscape = network.generate_cytoscape_json()
    assert "format_version" in cytoscape
    assert "elements" in cytoscape
    assert "nodes" in cytoscape["elements"]
    assert "edges" in cytoscape["elements"]

def test_save_cytoscape():
    dds_files_path = 'dummy_data/dummy_dds.yml'
    network = cn.new_network(dds_files_path,
                             save_name='dummy_results/network_test.cyjs')
    network.collect_genes()
    network.get_network()
    network.generate_cytoscape_json()
    network.save_network()
    # assert that the file save_name exists
    assert os.path.exists('dummy_results/network_test.cyjs')