### This file is the mail pipeline for my priject
import pandas as pd
import pytest
import yaml
import Constants
import cytoscape as ct
import network_processing as ntp
import os.path
from os import path
import py4cytoscape as py4
import node_evaluation as ne
import argparse



def get_cytoscape_network(file_name):
    cytoscape_json = ct.read_cytoscape_json(cytoscape_file=file_name)
    nodes, edges_metadata, relationship = ct.format_cytoscape_json(cytoscape_json=cytoscape_json)
    the_network = ntp.create_graph_from_dictionaries(nodes=nodes, edges=edges_metadata, relationship=relationship)
    return the_network


def add_mirnas_n_select(cytoscape_network, name, use_prefix=True, add_mirnas=True, cutoff=0.5):
    """

    :param cytoscape_network: str Name of the network saved as cyjs
    :param name: str Name to save the graph and subproducts
    :param use_prefix: Bool
    :param add_mirnas: Bool
    :return:
    """
    if use_prefix:
        px1 = "graph0_"
        px2 = "graph1_"
    else:
        px1 = px2 = ""
    if not os.path.exists(f"{px1}{name}.pkl"):
        the_network = get_cytoscape_network(cytoscape_network)
        ntp.save_graph(the_network, f"{px1}{name}.pkl")
    else:
        the_network = ntp.load_graph(f"{px1}{name}.pkl")
    if add_mirnas:
        if not os.path.exists(f"{px2}{name}.pkl"):
            ntp.add_mirna_relationships(the_network)
            ntp.save_graph(the_network, f"{px2}{name}.pkl")
        else:
            the_network = ntp.load_graph(f"{px2}{name}.pkl")

    network = ntp.remove_nodes_low_centrality(graph=the_network, cutoff=cutoff)
    ntp.set_positions(network)
    ntp.save_graph(network, f"{px2}_filtered_{name}.pkl")
    save_as_cjsn(network, f'{px2}{name}.cyjs')
    return network


def add_mirnas_n_tissues(cytoscape_network, name, use_prefix=True, add_mirnas=True, add_tissues=True, add_system=True):
    """
    :param add_system: Bool
    :param cytoscape_network: str Name of the network saved as cyjs
    :param name: str Name to save the graph and subproducts
    :param use_prefix: Bool
    :param add_mirnas: Bool
    :param add_tissues: Bool
    :return:
    """
    if use_prefix:
        px1 = "graph0_"
        px2 = "graph1_"
    else:
        px1 = px2 = ""
    the_network = get_network(px1, px2, cytoscape_network)
    if add_tissues:
        ntp.add_tissue_relationship(the_network)
    if add_system:
        ntp.add_organ_system_relationship(the_network)
    the_network.remove_nodes_from(['Tnf'])
    ntp.set_positions(the_network)
    ne.evaluate_nodes(the_network)
    # ne.remove_nodes(the_network, threshold=0.85)
    ntp.save_graph(the_network, f"{px2}{name}.pkl")
    save_as_cjsn(the_network, f'{px2}{name}.cyjs')


def save_as_cjsn(network, name):
    """
    This function will take a network and save it for cytoscape use
    :return:
    """
    json_network = ntp.convert_to_json(network)
    ct.save_cytoscape_json(json_network, cytoscape_file_name=name)


def get_network(px1, px2, cytoscape_network, save_name, add_mirnas=True):
    """
    This function will get the network, first is going to check if it already
    constructed it or if it has to calculate it.
    :param px1: the prefix of the first buiid
    :param px2: the prefix of the network with the microRNAS
    :param cytoscape_network: The name of the file that holds the cjsn of the network
    :param save_name: The name from which the network is saved
    :param add_mirnas: If we shoul add mirnas or not
    :return:
    """
    if not os.path.exists(f"{px1}{save_name}.pkl"):
        the_network = get_cytoscape_network(cytoscape_network)
        ntp.save_graph(the_network, f"{px1}{save_name}.pkl")
    else:
        the_network = ntp.load_graph(f"{px1}{save_name}.pkl")
    if add_mirnas:
        if not os.path.exists(f"{px2}{save_name}.pkl"):
            ntp.add_mirna_relationships(the_network)
            ntp.save_graph(the_network, f"{px2}{save_name}.pkl")
        else:
            the_network = ntp.load_graph(f"{px2}{save_name}.pkl")
    return the_network


def full_flow_pageRank(cytoscape_network, name, use_prefix=True, cutoff=0.5):
    """
    :param cytoscape_network: str Name of the network saved as cyjs
    :param name: str Name to save the graph and subproducts
    :param use_prefix: Bool
    :param cutoff:

    :return:
    """
    if use_prefix:
        px1 = "graph0_"
        px2 = "graph1_"
    else:
        px1 = px2 = ""
    the_network = get_network(px1, px2, cytoscape_network, save_name=name)

    network = ntp.remove_nodes_low_centrality_pageRank(graph=the_network, cutoff=cutoff)
    network = ntp.set_positions(network)
    ntp.save_graph(network, f"Networks_pkl/pageRank_{px2}_{name}.pkl")
    save_as_cjsn(network, f'Networks_CYJS(out)/pageRank_{px2}_{name}.cyjs')

    return network


def full_flow_genes_tf(cytoscape_network, name, use_prefix=True, dds_df=None, tissue_df=None, cutoff=0.5):
    """
    :param cytoscape_network: str Name of the network saved as cyjs
    :param name: str Name to save the graph and subproducts
    :param use_prefix: Bool
    :param cutoff:

    :return:
    """
    if use_prefix:
        px1 = "original_"
        px2 = "mirnas_"
    else:
        px1 = px2 = ""
    if not os.path.exists(f"Networks_pkl/metadata_{px2}_{name}.pkl"):
        network = get_network(px1, px2, cytoscape_network, save_name=name)
        ntp.add_pathways_to_nodes(graph=network,
                                  pathway_file="data/enr_pvals_RNAseq_abundances_adjusted_combat_inmose_young.vs"
                                               ".old_filtered.csv")
        ntp.mark_TF_nodes_from_file(graph=network, TF_file=)
        ntp.add_tissue_to_nodes(graph=network,
                                tissues_df=tissue_df)
        ntp.mark_miR_nodes(graph=network)
        if dds_df is not None:
            ntp.add_DDS_data(network, dds_df=dds_df)

        ntp.weight_nodes(graph=network)
        ntp.save_graph(network, f"Networks_pkl/metadata_{px2}_{name}.pkl")
    else:
        network = ntp.load_graph(f"Networks_pkl/metadata_{px2}_{name}.pkl")
    network = ntp.set_positions(network)

    network = ntp.remove_nodes_low_centrality_pageRank(graph=network, weigth='weigh', cutoff=cutoff)
    n_mirnas = ntp.get_n_mirs(network)
    #network = nx.get_interest_genes_and_neighbors(n_neighbors=2, graph= network)
    ntp.save_graph(network, f"Networks_pkl/complete_n_tf_{px2}_{name}.pkl")
    save_as_cjsn(network, f'Networks_CYJS(out)/complete_n_tf_{px2}_{name}.cyjs')

    return network


def open_cytoscape(file_name):
    py4.open_session(file_name)


def main(file, path_dds_data, path_tissue_data, ranks):
    if path_tissue_data is None:
        path_tissue_data = Constants.PATH_TISSUE_DATA
    if path_dds_data is None:
        path_dds_data = Constants.PATH_DDS_DATA
    if ranks is None:
        ranks = Constants.RANKS
    dds_df = pd.read_csv(path_dds_data, index_col=0).fillna(0)
    tissue_df = pd.read_csv(path_tissue_data, index_col=0).fillna(0)
    file_name = os.path.basename(file)
    
    for n in ranks:
        name = f"{file_name.split('.')[0]}_cutoff_{n}"
        full_flow_genes_tf(file, name, dds_df=dds_df, tissue_df=tissue_df, cutoff=n)

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process gene network analysis with different cutoff values.")
    parser.add_argument("--config", type=str, help="Path to the configuration YAML file.")
    parser.add_argument("--network", type=str, help="Path to the input network file.")
    parser.add_argument("--DE_data", type=str, help="Path to the normalized DeSeq2DataSet file.")
    parser.add_argument("--Tissue_expression_data", type=str, help="Path to the gene tissue data file.")
    parser.add_argument("--Cell_type_data", type=str, help="Path to the gene cell type data file.")
    parser.add_argument("--Pathway_enrichment_data", type=str, help="Path to the gene pathway enrichment data file.")
    parser.add_argument("--TF_enrichment_data",  type=float, help="Path to the gene TF enrichment data file.")
    parser.add_argument("--Cutoff", nargs='+', type=float, help="List of cutoff values for ranking separated by space.")
    parser.add_argument("--output", type=str, help="Path to the output directory.")
    
    args = parser.parse_args()
    if args.config:
        config_data = load_config(args.config)
    else:
        config_data = load_config("settings/metadata.yml")
    # Override config values with provided CLI arguments (if any)
    network = args.network if args.network else config_data.get("oNetwork")
    DE_data = args.DE_data if args.DE_data else config_data.get("path_DDS_data")
    Tissue_expression = args.Tissue_expression if args.Tissue_expression else config_data.get("path_tissue_data")
    Cell_type = args.Cell_type if args.Cell_type else config_data.get("path_cell_type_data")
    Pathway_enrichment = args.Pathway_enrichment if args.Pathway_enrichment else config_data.get("path_pathway_file")
    TF_enrichment = args.TF_enrichment if args.TF_enrichment else config_data.get("path_tf_file")
    Cutoff = args.Cutoff if args.Cutoff else config_data.get("page_rank_cutoff", [0.5, 0.7])
    print(f"Running main pipeline with the following parameters:\n"
          f"Network: {network}\n"
          f"DE_data: {DE_data}\n"
          f"Tissue_expression: {Tissue_expression}\n"
          f"Cell_type: {Cell_type}\n"
          f"Pathway_enrichment: {Pathway_enrichment}\n"
          f"TF_enrichment: {TF_enrichment}\n"
          f"Cutoff: {Cutoff}\n"
          f"Output: {args.output}")
    # main(config_data)

    