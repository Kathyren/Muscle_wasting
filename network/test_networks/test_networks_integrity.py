import sys
sys.path.append('../')
#from .. import network_processing as nx
#from .. import main_pipeline as mp
import pytest
import cytoscape as ct
from network import main_pipeline as mp
from network import network_processing
import yaml


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

