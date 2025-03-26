### This file is the mail pipeline for my project
import sys
import os.path

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pandas as pd
import pytest
import yaml
import cytoscape as ct
import network_processing as ntp

from os import path
import py4cytoscape as py4
import node_evaluation as ne
import argparse

# Set base directory as the MUSCLE_WASTING/


# Set the working directory to the main
os.chdir(BASE_DIR)




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


def full_flow_genes_tf(cytoscape_network, name, use_prefix=True, dds_df=None,
                       tissue_df=None, tf_file=None, pathway_file=None, cell_type=None,
                         coefficients={},cutoff=0.5, path = "network"):
    """
    :param pathway_file:
    :param tf_file:
    :param tissue_df:
    :param dds_df:
    :param cytoscape_network: str Name of the network saved as cyjs
    :param name: str Name to save the graph and subproducts
    :param use_prefix: Bool
    :param cutoff:

    :return:
    """
    # if any of the attributes is None, we are going to use the default values

    if use_prefix:
        px1 = "original_"
        px2 = "mirnas_"
    else:
        px1 = px2 = ""
    network=None
    if not os.path.exists(f"{path}Networks_pkl/metadata_{px2}_{name}.pkl"):

        network = get_network(px1, px2, cytoscape_network, save_name=name)
        #ntp.add_pathways_to_nodes(graph=network,
        #                          pathway_file=pathway_file)
        #ntp.mark_TF_nodes_from_file(graph=network, TF_file=tf_file)
        #ntp.add_tissue_to_nodes(graph=network,
        #                        tissues_df=tissue_df)
        ntp.mark_miR_nodes(graph=network)
        if dds_df is not None:
            dds_df = dds_df.div(dds_df.sum(axis=1), axis=0)
            ntp.add_DDS_data(network, dds_df=dds_df)
        if tf_file is not None:
            #ntp.mark_TF_nodes_from_file(graph=network, TF_file=tf_file) 0-1
            tf_df = pd.read_csv(tf_file, index_col=0)
            # normalize the dataframe

            #tf_df = (tf_df - tf_df.min()) / (tf_df.max() - tf_df.min())
            ntp.add_TF_data_from_df(graph=network,tf_df=tf_df)
        if tissue_df is not None:
            ntp.add_other_data(graph=network, data_df=tissue_df, name='tissue')
        if pathway_file is not None:
            # normalize the dataframe
            pathway_df = pd.read_csv(pathway_file, index_col=0)
            ntp.add_pathways_to_nodes_from_df(graph=network, all_pathway_df=pathway_df)
        if cell_type is not None:
            
            ntp.add_other_data(graph=network, data_df=cell_type, name='cell_type')

        print(f"Saving on {path}Networks_pkl/metadata_{px2}_{name}.pkl  reading it")
        ntp.save_graph(network, f"{path}Networks_pkl/metadata_{px2}_{name}.pkl")
    else:
        print(f"The nerwork is saved on {path}Networks_pkl/metadata_{px2}_{name}.pkl  reading it")
        network = ntp.load_graph(f"{path}Networks_pkl/metadata_{px2}_{name}.pkl")

    #if dds_df is not None:
    #        ntp.add_DDS_data(network, dds_df=dds_df)
    #if tf_file is not None:
    #    #ntp.mark_TF_nodes_from_file(graph=network, TF_file=tf_file)
    #    ntp.add_TF_data_from_file(graph=network,tf_file=tf_file)
    #if tissue_df is not None:
    #    ntp.add_other_data(graph=network, data_df=tissue_df, name='tissue')
    #if pathway_file is not None:
    #    ntp.add_pathways_to_nodes(graph=network, pathway_file=pathway_file)
    #if cell_type is not None:
    #    ntp.add_other_data(graph=network, data_df=cell_type, name='cell_type')
    network = ntp.set_positions(network)
    ntp.weight_nodes(graph=network, coefficients=coefficients)
    if 'miR_enhancement' in coefficients.keys():
        ntp.weight_edges(graph=network, node_weight='weigh',
                         mir_enhancer= coefficients['miR_enhancement'])
    else:
        ntp.weight_edges(graph=network, node_weight='weigh')

    network = ntp.remove_nodes_low_centrality_pageRank(graph=network, weight='weightScore', cutoff=cutoff)
    #n_mirnas = ntp.get_n_mirs(network)
    #network = nx.get_interest_genes_and_neighbors(n_neighbors=2, graph= network)
    ntp.save_graph(network, f"network/Networks_pkl/complete_n_tf_{px2}_{name}.pkl")
    ntp.save_graph(network, f"{path}/{name}.pkl")
    ntp.separate_metadata(graph=network )
    cytoscape_file = f'{path}Networks_CYJS(out)/complete_n_tf_{px2}_{name}.cyjs'
    save_as_cjsn(network, cytoscape_file)
    cytoscape_file = f'{path}/{name}.cyjs'
    save_as_cjsn(network, cytoscape_file)

    return network, cytoscape_file


def open_cytoscape(file_name):
    py4.open_session(file_name)


def main(config_data, file, path_dds_data, path_tissue_data, pathway_file,
          tf_file, cell_type_file, ranks, open_cytoscape_=False, name_output=None, path= "network/"):
    
    if file is None:
        file = config_data.oNetwork
    if path_tissue_data is None:
        path_tissue_data = config_data.path_tissue_data
    if path_dds_data is None:
        path_dds_data = config_data.path_DDS_data
    if pathway_file is None:
        pathway_file = config_data.path_pathway_file
    if tf_file is None:
        tf_file = config_data.path_tf_file
    if ranks is None:
        ranks = config_data.page_rank_cutoff
    if cell_type_file is None:
        cell_type_file = config_data.path_cell_type_data

    dds_df = pd.read_csv(path_dds_data, index_col=0).fillna(0)
    tissue_df = pd.read_csv(path_tissue_data, index_col=0).fillna(0)
    print(f"cell_type_file: {cell_type_file}")
    #print(f"config_data.path_cell_type_data: {config_data.path_cell_type_data}")

    cell_type_df = pd.read_csv(cell_type_file, encoding="latin1", index_col=0).fillna(0)
    if name_output is None:
        file_name = os.path.basename(file)
    else:
        file_name = name_output
    
    for n in ranks:
        name = f"{file_name.split('.')[0]}_cutoff_{n}"
        network, cytoscape_file = full_flow_genes_tf(file, name, 
                           dds_df=dds_df, 
                           tissue_df=tissue_df,
                           tf_file=tf_file, 
                           pathway_file=pathway_file,
                           cell_type= cell_type_df,
                           coefficients=config_data.get('coefficients'),
                           cutoff=n, path=path)
        if open_cytoscape_:
            open_cytoscape(cytoscape_file)
    return network

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
    parser.add_argument("--Pathway_enrichment", type=str, help="Path to the gene pathway enrichment data file.")
    parser.add_argument("--TF_enrichment",  type=float, help="Path to the gene TF enrichment data file.")
    parser.add_argument("--Cutoff", nargs='+', type=float, help="List of cutoff values for ranking separated by space.")
    parser.add_argument("--network_name", type=str, help="Name of the network.")
    parser.add_argument("--output", type=str, help="Path to the output directory.")
    parser.add_argument("--open_cytoscape", type=bool, help="Open the Cytoscape session after processing.", default=False)
    
    args = parser.parse_args()
    if args.config:
        config_data = load_config(args.config)
    else:
        config_data = load_config("network/settings/metadata.yml")
    # Override config values with provided CLI arguments (if any)
    network = args.network if args.network else config_data.get("oNetwork")
    DE_data = args.DE_data if args.DE_data else config_data.get("path_DDS_data")
    Tissue_expression = args.Tissue_expression_data if args.Tissue_expression_data else config_data.get("path_tissue_data")
    Cell_type = args.Cell_type_data if args.Cell_type_data else config_data.get("path_cell_type_data")
    Pathway_enrichment = args.Pathway_enrichment if args.Pathway_enrichment else config_data.get("path_pathway_file")
    TF_enrichment = args.TF_enrichment if args.TF_enrichment else config_data.get("path_tf_file")
    Cutoff = args.Cutoff if args.Cutoff else config_data.get("page_rank_cutoff", [0.5, 0.7])
    output = args.output if args.output else config_data.get("output")
    open_cytoscape_ = args.open_cytoscape if args.open_cytoscape else config_data.get("open_cytoscape", False)
    coef = config_data.get("coefficients")
    name = args.network_name if args.network_name else os.path.basename(network)
    print(coef)
    print(f"Running main pipeline with the following parameters:\n"
          f"Network: {network}\n"
          f"DE_data: {DE_data}\n"
          f"Tissue_expression: {Tissue_expression}\n"
          f"Cell_type: {Cell_type}\n"
          f"Pathway_enrichment: {Pathway_enrichment}\n"
          f"TF_enrichment: {TF_enrichment}\n"
          f"Cutoff: {Cutoff}\n"
          f"Output: {output}\n"
          f"Open Cytoscape: {open_cytoscape_}")
    main(config_data=config_data, file=network, path_dds_data=DE_data,
          path_tissue_data=Tissue_expression, pathway_file=Pathway_enrichment, 
          tf_file=TF_enrichment, cell_type_file=Cell_type, ranks=Cutoff,
         open_cytoscape_=open_cytoscape_, name_output=name)


    