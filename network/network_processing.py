import pwd
import random

import networkx
import networkx as nx
import numpy as np
import pandas as pd
from networkx import neighbors
import sys
import os

from pandas.core.interchange.dataframe_protocol import DataFrame

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from database_analysis import sql_operations as sql
from cytoscape import read_cytoscape_json, format_cytoscape_json, create_cytoscape_node, create_cytoscape_edge, protein_name


def create_graph(mirnas, genes, relationsip):
    """
    This function will receive a list of genes and mirnas and create the network of them
    :param mirnas: List of string with the mirna names
    :param genes: List of string with the gene names
    :param relationsip: List of tuple of string with the node names
    """
    G = nx.Graph()
    G.add_nodes_from(mirnas, size=10, type="mirna")
    G.add_nodes_from(genes, type="gene")
    G.add_edges_from(relationsip)
    is_bipartite(G)
    return G


def create_graph_from_dictionaries(nodes, relationship, edges):
    """
    This function will receive a list of genes and mirnas and create the network of them
    :param edges: List of dictionaries with the metadata of each relationship
    :param mirnas: List of dictionary with the mirna names
    :param genes: List of string with the gene names
    :param relationsip: List of tuple of string with the node names
    """

    #protein_name = source.get_main_name()
    G = nx.DiGraph()

    for element in nodes:
        node = element['data'][protein_name]
        if node not in G._node:
            G.add_node(node, **element)
    for edge in edges:
        try:
            if edge['target'] not in G._adj[edge['source']]:
                G.add_edge(edge['source'], edge['target'], **edge)
        except Exception as e:
            print(e)

    # is_bipartite(G)

    return G


def draw_graph(G):
    """
    This function will draw the graph in a matplotlib window. It is not working... I think
    :param G: The networkx object
    :return: None
    """
    color_map = []

    for node, data in G.nodes(data=True):
        if data['type'] == 'mirna':
            # color map pink
            #color_map.append()
            color_map.append(0.25)  # blue color
        elif data['type'] == 'gene':
            color_map.append(0.7)  # yellow color
    nx.draw(G, node_color=color_map, with_labels=False)
    # nx.draw(G)


def is_bipartite(G):
    return True



def filter_by_degree(G):
    """
    This function will filter the nodes by degree
    :param G: The networkx object
    :return: The networkx object with the nodes filtered
    """
    graph = G.copy()
    centrality = nx.degree(G)
    c = list(centrality)
    _, values = zip(*c)
    average = np.mean(values)
    sandar_dev = np.std(values)
    x = np.quantile(values, 0.75)
    delete_nodes = []
    for node in centrality:
        if node[1] < x:
            delete_nodes.append(node[0])
    graph.remove_nodes_from(delete_nodes)
    save_graph(graph, "high_degree.pkl")
    return graph


def remove_nodes_low_centrality(graph, cutoff=0.75):
    """
    This function will calculate a centrality


    :param cutoff: What is going to be the quartile of nodes that we are going to keep
    :param graph: A networkX graphThe graph without the unrelevant nodes
    :return:
    """
    # @Todo: Establish a centrality or make this dynamic in code
    # dc = list(nx.degree_centrality(graph).values())
    # idc = nx.in_degree_centrality(G)
    # odc = nx.out_degree_centrality(G)
    # bc = list(nx.betweenness_centrality(graph).values())
    # lc = list(nx.load_centrality(graph).values())
    #    ec = list(nx.eigenvector_centrality(G).values())
    cc_t = nx.closeness_centrality(graph)
    cc = list(cc_t.values())
    x = np.quantile(cc, cutoff)
    delete_nodes = []
    delete_edges = []
    for node in cc_t.items():
        if node[1] < x:
            delete_nodes.append(node[0])
            delete_edges.extend(graph.edges(node[0]))
    graph.remove_edges_from(delete_edges)
    graph.remove_nodes_from(delete_nodes)
    # draw_graph(graph)

    return graph


def evaluate_edge_quality(graph):
    """
    This function will evaluate if the edges have all the necesary requierement

    :param graph:
    :return:
    """
    #todo: Change this dictionary with a ymal file eventually
    required_keys = {"id", "interaction", "name", "selected", "shared_interaction", "shared_name", "source",
                     "target", "weight"}
    nodes_with_ids = {data['data']['id']: node for node, data in graph.nodes(data=True) if
                      'data' in data and isinstance(data['data'], dict) and 'id' in data['data']}
    reverse_dict = {v: k for k, v in nodes_with_ids.items()}

    edges_to_remove = []

    for u, v, data in graph.edges(data=True):
        if 'data' in data:
            if isinstance(data['data'], dict):
                missing_keys = required_keys - data['data'].keys()
                for key in missing_keys:
                    data['data'][key] = None  # or any default value you prefer
                source_id = data['data'].get('source')
                target_id = data['data'].get('target')

                if source_id not in nodes_with_ids or target_id not in nodes_with_ids:
                    remove = False
                    try:
                        if source_id not in nodes_with_ids:
                            if u in nodes_with_ids.values():
                                source_id = reverse_dict.get(u)
                                data['data']['source'] = source_id
                            else:
                                remove = True
                        if target_id not in nodes_with_ids:
                            if v in nodes_with_ids.values():
                                target_id = reverse_dict.get(v)
                                data['data']['target'] = target_id
                            else:
                                remove = True
                        if remove:
                            edges_to_remove.append((u, v))

                    except Exception as e:
                        edges_to_remove.append((u, v))

            else:
                edges_to_remove.append((u, v))
        else:
            edges_to_remove.append((u, v))

    graph.remove_edges_from(edges_to_remove)


def evaluate_node_quality(graph):
    """
    This function will evaluate if the nodes have all the necesary requierement

    :param graph:
    :return:
    """
    #todo: Change this dictionary with a ymal file eventually
    required_keys = {"id", "SUID", "name", "selected", "shared_name"}

    nodes_to_remove = []

    for node, data in graph.nodes(data=True):
        if 'data' in data:
            if isinstance(data['data'], dict):
                missing_keys = required_keys - data['data'].keys()
                for key in missing_keys:
                    data['data'][key] = None  # or any default value you prefer
            else:
                nodes_to_remove.append(node)
        else:
            nodes_to_remove.append(node)

    graph.remove_nodes_from(nodes_to_remove)


def remove_node_and_edges(graph, node):
    """
        This function will calculate a centrality
        :param cutoff: What is going to be the quartile of nodes that we are going to keep
        :param graph: A networkX graphThe graph without the unrelevant nodes
        :return:
    """

    delete_nodes = [node]
    delete_edges = []
    #delete_edges.extend(graph.edges(node))
    delete_edges.extend(graph.in_edges(node))
    delete_edges.extend(graph.out_edges(node))
    graph.remove_edges_from(delete_edges)
    #graph.remove_nodes_from(delete_nodes)
    # draw_graph(graph)

    return graph


def mark_TF_nodes_from_file(graph, TF_file):
    """
    This function takes a fiile that cointains the TF and a score value to tell us how relevant they are
    TF, enrichment and pvals

    This will get the threshold feature, take the TF that pass the filter and give a dictionary with the nodes as
    key and the TF as values

    :param pathway_file:
    :return:
    """
    if TF_file is None:
        return None
    tf_dic = pd.read_csv(TF_file, index_col=0).to_dict('index')

    for node, data in graph.nodes(data=True):
        if 'data' in data and 'id' in data['data'] and data['data']['name'] in tf_dic.keys():
            graph.nodes[node]['data']['node_type'] = 'TF'
            graph.nodes[node]['data']['presence'] = tf_dic[data['data']['name']]['enrichment']
            graph.nodes[node]['data']['TF_enrrichment'] = tf_dic[data['data']['name']]['enrichment']

    pass

def add_TF_data_from_file(graph, tf_file):
    """
    This function takes a fiile that cointains the TF and a score value from each condition and 
    will add them to their respective nodes, similar to add_DDS_data
    It will add in metadata the dictionary of conditions and the values of the TF

    :param graph:
    :param TF_file: The file with the TF and the score values. The index is the TF name and the 
    columns are the conditions. 
    :return:
    """
    if tf_file is None:
        return None
    tf_df = pd.read_csv(tf_file, index_col=0)
    for index, row in tf_df.iterrows():
        add_tf_to_node(graph=graph, node_name=index, tf=dict(row))

    pass

def add_tf_to_node(graph, node_name, tf: dict):
    """
    This function will add as metadata to the node a list of tf
    :param graph:
    :param node:
    :param dds: A dictionary with the Enriched score TF values
    :return:
    """
    for node, data in graph.nodes(data=True):
        if 'data' in data and 'id' in data['data'] and data['data']['name'] == node_name:
            for tf, value in tf.items():
                if not 'metadata' in graph.nodes[node]['data']:
                        graph.nodes[node]['data']['metadata']={}
                if 'tf' in graph.nodes[node]['data']['metadata']:
                    graph.nodes[node]['data']['metadata']['tf'][tf]=value
                else:
                    graph.nodes[node]['data']['metadata']['tf'] = {tf:value}
                

def mark_miR_nodes(graph):
    """
    This function takes a fiile that cointains the TF and a score value to tell us how relevant they are
    TF, enrichment and pvals

    This will get the threshold feature, take the TF that pass the filter and give a dictionary with the nodes as
    key and the TF as values

    :param pathway_file:
    :return:
    """
    for node, data in graph.nodes(data=True):
        if str('data' in data and 'id' in data['data'] and data['data']['name']).startswith('hsa-miR'):
            graph.nodes[node]['data']['node_type'] = 'miR'
            graph.nodes[node]['data']['presence'] = None

    pass


def get_n_mirs(graph):
    """
    This function takes a fiile that cointains the TF and a score value to tell us how relevant they are
    TF, enrichment and pvals

    This will get the threshold feature, take the TF that pass the filter and give a dictionary with the nodes as
    key and the TF as values

    :param pathway_file:
    :return:
    """
    mir_nodes = 0
    for node, data in graph.nodes(data=True):
        if 'data' in data and 'node_type' in data['data'] and data['data']['node_type'] == 'miR':
            mir_nodes += 1

    return mir_nodes

def get_mirs(graph):
    """
    This function takes a fiile that cointains the TF and a score value to tell us how relevant they are
    TF, enrichment and pvals

    This will get the threshold feature, take the TF that pass the filter and give a dictionary with the nodes as
    key and the TF as values

    :param pathway_file:
    :return:
    """
    mir_nodes = []
    for node, data in graph.nodes(data=True):
        if 'data' in data and 'node_type' in data['data'] and data['data']['node_type'] == 'miR':
            mir_nodes.append(node)

    return mir_nodes

def get_target_of_mir(G, mir):
    """
    Goes to the network and gets the targets of a particular mir
    :param G:
    :param mir:
    :return:
    """
    return list(G.successors(mir))

def add_DDS_data(graph, dds_df):
    """
    This function will look for the nodes in the gene column and then add as an attribute every of the columns
    :param graph:
    :param dds_df:
    :return:
    """
    for index, row in dds_df.iterrows():
        add_dds_to_node(graph=graph, node_name=index, dds=dict(row))

    pass


def add_dds_to_node(graph, node_name, dds: dict):
    """
    This function will add as metadata to the node a list of pathways
    :param graph:
    :param node:
    :param dds: A dictionary with the Dseq2DataSet values (stat or log2FoldChange)
    :return:
    """
    for node, data in graph.nodes(data=True):
        if 'data' in data and 'id' in data['data'] and data['data']['name'] == node_name:
            for dd, value in dds.items():
                if not 'metadata' in graph.nodes[node]['data']:
                        graph.nodes[node]['data']['metadata']={}
                if 'dds' in graph.nodes[node]['data']['metadata']:
                    graph.nodes[node]['data']['metadata']['dds'][dd]=value
                else:
                    graph.nodes[node]['data']['metadata']['dds'] = {dd:value}
                graph.nodes[node]['data'][dd] = value


def add_metadata_to_node(graph, node_name, tissues: dict, name='tissue_expr'):
    """
    This function will add the tissue data to the node.
    :param graph: The networkx graph object
    :param node_name: The name of the node to get the tissue. This is the name of the gene
    :param tissues: A dictionary with the tissues and the expression. The key should be the gene as in the node.
    :param name: The name of the attribute, tissue is the default (tissue_expr)
    :return:
    """
    for node, data in graph.nodes(data=True):
        if 'data' in data and 'id' in data['data'] and data['data']['name'] == node_name:
            if not 'metadata' in graph.nodes[node]['data']:
                graph.nodes[node]['data']['metadata'] = {}
            graph.nodes[node]['data']['metadata'][name] = tissues


def add_tissue_to_nodes(graph, tissues_df):
    """
    This function will get the data in tissues_df and add them to the node
    Is expected that tissues_df is only the tissues, and it is already processed
    (the normalize Transcripts per million of the tissue divided by the sum of the normalized Transcripts per Million)
    :param graph:
    :param tissues_df:
    :return:
    """
    filtered_tissues_df = tissues_df[tissues_df.index.isin(graph.nodes)]
    for index, row in filtered_tissues_df.iterrows():
        add_metadata_to_node(graph=graph, node_name=index, tissues=dict(row))


def extract_genes_from_pathways(pathway_df:DataFrame, threshold_feature='Combined score',
                                id_feature='Features', pathway_feature='Term', threshold_value=50):
    """
    Soon to be deprecated
    This function takes a fiile that cointains
    id term set size overlap ratio p_value fdr p-vale odds ratio combine score and features


    :param threshold_feature:
    :param id_feature:
    :param pathway_feature:
    :param pathway_file:
    :param threshold_value:
    :return:
    """
    #pathway_df = pd.read_csv(pathway_file)
    #pathway_df = pathway_df[pathway_df[threshold_feature] > threshold_value]
    feature_dict = {}

    # Iterate through each row in the dataframe
    for index, row in pathway_df.iterrows():
        # Split the Features column by ';' to get individual features
        features = row[id_feature].split(',')
        features = row[id_feature].replace("'", "").replace(" ","")
        # Now features is a string like "['a','b','c']". COnvert to the list
        features = features[1:-1].split(',')
        # Iterate through each feature
        for feature in features:
            if feature:  # Ignore empty strings
                if feature not in feature_dict:
                    feature_dict[feature] = []
                # Append the Term (pathway) to the list of pathways for this feature
                feature_dict[feature].append(row[pathway_feature])
    return feature_dict


def add_pathways_to_nodes(graph, pathway_file):
    """
    This function will add the pathway to the node. Currently it takes the string directory where the pathways are.
    It will call the funtion extract_genes_from_pathways to get the pathways and then add them to the nodes to 
    generate the dictorionary of the pathways. That dictionary is added to the node with add_pathway_to_node

    :param graph: The networkx graph object with the format of mirKat network (see Documentation)
    :param pathway_file: Path to the file with the pathways. See Documentation for the format
    :return: None

    """
    all_pathway_df = pd.read_csv(pathway_file)
    gene_pathway_dic = extract_genes_from_pathways(all_pathway_df, threshold_feature='Combined score',
                                                   id_feature='genes', pathway_feature='Term', threshold_value=10)
    for gene, pathways in gene_pathway_dic.items():
        _, svd, _ = get_SVD_pathways(pathway_df=all_pathway_df.drop(columns=['genes']),
                                     pathways_list=pathways)
        add_pathway_to_node(graph, gene, pathways, svd=svd)
    return graph

def get_SVD_pathways(pathway_df:DataFrame, pathways_list:list, pathway_column = 'Term'):
    """
     # Pathways have the list of pathways that the gene belongs
        # the pathway file will have the dataframe of the combinations
        # and the combined score of the pathways on each combination.
        # I can use pathways to get the sub dataframe of the pathway_file
        # and use that to calculate the Singular Value Decomposition (SVD)
        # and add that value

    This function will get the SVD of the pathways in the list.
    Will get the subset of pathways form pathway list (holded in "term" column) 
    and calculate the SVD of the subset of pathways.

    :param pathway_df: The dataframe with the pathways. The genes column should be already removed
    :param pathways_list: The list of pathways to calculate the SVD
    :param pathway_column: The column where the pathways are
    :return: The singular vector of the pathways,
      the singular value and the singular vector of the conditions
    """
    subset_pathways = pathway_df[pathway_df[pathway_column].isin(pathways_list)]
    subset_pathways.set_index(pathway_column, inplace=True)
    U, Sigma, VT = np.linalg.svd(subset_pathways, full_matrices=False)
    # Sigma is returned as a 1D array; convert to diagonal matrix
    # Sigma[0] the largest singular value (dominant pattern)
    # U[:,0] First pathway singular vector (importance of each pathway)
    # VT[0,:] First condition singular vector (importance of each condition)
    return U[:,0], Sigma[0], VT[0,:]



def add_pathway_to_node(graph, node_name, pathways, svd = 0):
    """
    This function will add as metadata to the node a list of pathways
    :param graph:
    :param node:
    :param pathways:
    :return:
    """
    for node, data in graph.nodes(data=True):
        if ('data' in data and 'id' in data['data'] and
                data['data']['name'] == node_name):
            graph.nodes[node]['data']['pathways'] = pathways
            if not 'metadata' in graph.nodes[node]['data']:
                graph.nodes[node]['data']['metadata'] = {}
            graph.nodes[node]['data']['metadata']['pathways'] = pathways
            graph.nodes[node]['data']['metadata']['pathways_svd'] = svd

def get_interest_genes_and_neighbors(graph, n_neighbors: int, distance: int, interest_genes: list):
    """
    This function will check a list of genes of interest and only keep those genes, and their nodes at distance n_neighbors
    :param graph:
    :param n_neighbors:
    :return:
    """

    for gene in interest_genes:
        pass

    pass


def select_next_valid_node(i, path_length, graph, neighbors) -> str:
    """
    This function will select the next node to go in the random walk.
    The node must have at least one neighbor to continue the path.
    It will check that data property exists and that the node has 
    the pathways property.


    :param i:
    :param path_length:
    :param graph:
    :param neighbors:
    :return:
    """
    next_node = random.choice(neighbors)
    #return next_node
    if i < (path_length - 1):
        n_neighbors = len(list(graph.successors(next_node)))
        while len(neighbors) > 1 > n_neighbors:
            data = graph.nodes[next_node]
            if 'data' in data and 'pathways' in data['data']:
                break
            neighbors.remove(next_node)
            next_node = random.choice(neighbors)
            n_neighbors = len(list(graph.successors(next_node)))
    return next_node


def get_random_path(G, start_node, path_length):
    """
    This function will get a random path of a certain length from a starting node in a graph.
    The path will be a list of nodes in the order they were visited.
    
    :param G: The graph to be traversed. This node should have the format of mirKat network (see Documentation)
    :param start_node: The starting node of the path
    :param path_length: The length of the path

    :return: A list of nodes in the order they were visited including the start node

    """
    path = [start_node]
    current_node = start_node

    for i in range(path_length - 1):
        try:
            neighbors = list(G.successors(current_node))  # Get the successors (out-neighbors) of the current node

        except KeyError as ke:
            return []
        if not neighbors:
            break  # No more neighbors to explore, end the path here

        next_node = select_next_valid_node(i, path_length, G, neighbors)

        path.append(next_node)
        current_node = next_node

    return path


def random_walk(graph, node_name, distance, sample_size) -> list:
    """
    This function takes  the graph, looks for the correspongng node and do the random path.
    :param graph:
    :param node_name:
    :param distance:
    :param sample_size:
    :return:
    """
    if node_name is None:
        return []
    paths = []
    for i in range(sample_size):
        path = get_random_path(graph, node_name, distance)
        paths.append(path)
    return paths


def remove_nodes_low_centrality_pageRank(graph, weigth=None, cutoff=0.75):
    """
    This function will calculate a centrality


    :param cutoff: What is going to be the quartile of nodes that we are going to keep
    :param graph: A networkX graphThe graph without the un-relevant nodes
    :return:
    """
    # @Todo: Establish a centrality or make this dynamic in code
    # dc = list(nx.degree_centrality(graph).values())
    # idc = nx.in_degree_centrality(G)
    # odc = nx.out_degree_centrality(G)
    # bc = list(nx.betweenness_centrality(graph).values())
    # lc = list(nx.load_centrality(graph).values())
    #    ec = list(nx.eigenvector_centrality(G).values())
    cc_t = nx.pagerank(graph, weight=weigth)
    cc = list(cc_t.values())
    x = np.quantile(cc, cutoff)
    delete_nodes = []
    delete_edges = []
    for node in cc_t.items():
        if node[1] < x:
            delete_nodes.append(node[0])
            delete_edges.extend(graph.edges(node[0]))
    graph.remove_edges_from(delete_edges)
    graph.remove_nodes_from(delete_nodes)
    # draw_graph(graph)

    return graph


def calculate_data(centralities):
    """
    This function takes a bounch of data (orginiginally intended as diferent centralities) and calculates how each
    centrality change among nodes to originally help to decide wich centrality to use
    :param centralities: = list of centralities

    :return: A dictionary of the centralities average, std, max and min value
    """

    av = []
    std = []
    max = []
    min = []
    for c in centralities:
        av.append(np.mean(c))
        std.append(np.std(c))
        max.append((np.max(c)))
        min.append(np.min(c))

    metrics = {'average': av,
               'std': std,
               'max': max,
               'min': min}
    return metrics


def save_graph(G, name):
    """
    This function will save a network obj to pickle file
    :param G:
    :param name:
    :return:
    """
    if int(nx.__version__[0]) > 2:
        import pickle
        with open(name, 'wb') as f:
            pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    else:
        nx.write_gpickle(G, name)


def convert_to_json(G):
    """
    This function will take a networkx object to json
    :param G:
    :return: json object with
    """
    json = nx.cytoscape_data(G)
    json = nx.node_link_data(G)
    json['elements'] = {'nodes': json.pop('nodes'), 'edges': json.pop('links')}

    return json


def load_graph(name):
    """
    Get a graph from a pickle file
    :param name:
    :return:
    """
    if int(nx.__version__[0]) > 2:
        import pickle
        with open(name, 'rb') as f:
            G = pickle.load(f)
    else:
        G = nx.read_gpickle(name)
    return G


get_random_relationships = "Select Distinct mrna, mirna_mature from binding limit 10;"
get_muscle_relationships = 'Select Distinct mrna, mirna_mature from binding where  ' \
                           'mirna_mature like "hsa-%" and mrna in (select gene from gene_bank) ;'


def get_nodes_names(G):
    """
        Return the name of the nodes in a list of string
    :param G:
    :return:
    """
    nodes = G.nodes
    return nodes


def get_mirna_mrna_relationships(genes):
    """
    This function will take the gene names and get the list
    of microrna and their relationships

    :param genes: list of string; List of gene names
    :return: The list of tuples for the relationships and the list of string with the name of the mirnas and the scores
    """
    genes = '"' + '","'.join(genes) + '"'
    query = 'Select Distinct mrna, mirna_mature from gene_mirna where  ' \
            f'mirna_mature like "hsa-%" and mrna in ({genes});'
    relationships = sql.get_query(query=query)
    genes = list(relationships['mrna'])
    mirnas = list(relationships['mirna_mature'])
    relationship = list(zip(mirnas, genes))
    scores = [1] * len(genes)  # list(relationships['probability'])

    return relationship, mirnas, scores


def get_mirna_tissue_edges(mirnas):
    """
    This function will collect the organ-mirna edges

    :return: tuple of the relationships and their probability
    """
    mirna_list = '"' + '","'.join(mirnas) + '"'
    query = f"select * from mirna_tissues where auto_mature in ({mirna_list});"
    result = sql.get_query(query=query)
    organ = list(result['organ'])
    mirnas = list(result['auto_mature'])
    relationship = list(zip(organ, mirnas))
    scores = list(result['organ_TSI'])
    return relationship, organ, scores


def get_tissue_system_edges(organs):
    """
        This function will collect the organ-system edges

        :return: tuple of the relationships and the systems
        """
    organ_list = '"' + '","'.join(organs) + '"'
    query = f"select distinct mirna_tissues.system, organ from mirna_tissues " \
            f"where organ in ({organ_list});"
    result = sql.get_query(query=query)
    organ = list(result['organ'])
    systems = list(result['system'])
    relationship = list(zip(organ, systems))
    return relationship, systems


def create_my_graph(query):
    """
    This function will create a graph from a query in the database
    :param query: The query to get the relationships
    :return: The networkx object
    """

    relationships = sql.get_query(query=query)
    genes = list(relationships['mrna'])
    mirnas = list(relationships['mirna_mature'])
    relationship = list(zip(genes, mirnas))
    graph = create_graph(mirnas=mirnas, genes=genes, relationsip=relationship)
    return graph


def add_relationships(graph, relationships):
    """
    This function will add to an existen network the relationships
    :param graph: The networkx objs
    :param relationships: A list of tuples with the relatioships
    :return: None
    """
    graph.add_edges_from(relationships)


def set_positions(network):
    """
    This function will rearrange the nodes and update the position as specified by cytoscape
    :param network:
    :return:
    """
    evaluate_node_quality(network)
    network_copy = network.copy()
    pos = nx.spectral_layout(network, scale=100, center=[50, 50])
    # pos = nx.circular_layout(network, scale=100, center=[50, 50])

    for node in network.nodes:
        try:
            network_copy.nodes[node]['position']['x'] = pos[node][0]
            network_copy.nodes[node]['position']['y'] = pos[node][1]

        except Exception as e:
            print(f"error with node {node} with {str(e)}")

            #network_copy.nodes[node]['position'] = {'x':pos[node][0],'y':f'{pos[node][1]}' }
            #network.nodes[node]['data'] = {'id': '', 'share_name': 'name', 'selected': ''}
            network_copy = remove_node_and_edges(network_copy, node)

            #network.node_evaluation.remove_nodes(node)
        #if not 'data' in network_copy.nodes[node]:
        #network_copy.nodes[node]['data'] = {'id': node, 'share_name': node, 'selected': False}

    evaluate_edge_quality(network_copy)

    #for edge in network_copy.edges:

    #    network_copy.edges[edge]['data'] = {'id': '', 'share_name': 'name', 'selected': ''}

    return network_copy


def add_mirna_relationships(network):
    """
    This function will add the mirna relationships to the network.
    It will get the node names, that at the point are the gene names, and get the mirna relationships
    through a sql call done in get_mirna_mrna_relationships. Up to now, the query is stablished and \
    cannot be changed but it is expected to be changed in the future.
    It will use the relationships to add the nodes (mirnas) and the edges to the network. The weigths of the edges
    are set to -1. 
    
    :param network: The networkx object
    :return: None
    
    """


    genes = get_nodes_names(network)
    relationships, mirnas, scores = get_mirna_mrna_relationships(genes)
    for idx, mirna in enumerate(mirnas):
        node_data = create_cytoscape_node(node_name=mirna, node_type='mirna', source='mirbase',
                                          node_data={'id': f'900{idx}', 'display_name': mirna, 'shared_name': mirna,
                                                     'name': mirna})
        network.add_node(mirna, **node_data)

    add_edge_from_relationships(network=network, edges=relationships, scores=scores)
    # add_relationships(network, relationships=relationships)


def add_edge_from_relationships(network, edges, scores=None, relationship_type="pm"):
    """
    This function will add the edges to the network.
     It will add the nodes if they are not present in the network. The edge properties are set to:
        shared_interaction, shared_name, interaction, selected, weight
    From those, weight is the most important since that will be used to calculate the influence of the nodes.



    :param network: The networkx object
    :param edges: The list of tuples with the relationships
    :param scores: The scores of the relationships
    :param relationship_type: The type of the relationship
    :return

    """

    # edges, scores = get_mirna_tissue_edges(mirnas=mirnas)
    if scores is None:
        scores = [1] * len(edges)
    for idx, edge in enumerate(edges):
        try:
            source = edge[0]
            target = edge[1]
            if 'data' not in network.nodes[source]:
                network.nodes[source]['data'] = {'id': source}
            if 'data' not in network.nodes[target]:
                network.nodes[target]['data'] = {'id': target}
            source_data = network.nodes[source]['data']
            target_data = network.nodes[target]['data']
            data2fill = {"id": f'800{idx}',
                         "source": source_data['id'],
                         "target": target_data['id'],
                         "shared_interaction": relationship_type,  ## predicted mirna
                         "shared_name": f"{source} ({relationship_type}) {target}",
                         "name": f"{source} ({relationship_type}) {target}",
                         "interaction": relationship_type,
                         "selected": False,
                         "weight": -1}
            edge_data = create_cytoscape_edge(source=source, target=target, node_data=data2fill, weight=scores[idx])
            network.add_edge(source, target, **edge_data)
        except Exception as e:
            print(f"Relationship {source}-{target} couldn't be added due to {str(e)}")
    add_relationships(network, edges)


def add_tissue_relationship(network):
    """
    This function will add the tissue relationships to the network.
    This function will take (for now) the tissues as a new node and add the edges to the network.
    The tissues are taken from the mirna_tissues table in the database. The relationships are added to the network
    with the add_edge_from_relationships function. The weigths of the edges are set to -1.

    This function will probably be retired. 


    :param network: The networkx object
    :return: None
    """

    nodes = get_nodes_names(network)
    relationships, organs, scores = get_mirna_tissue_edges(nodes)
    for idx, organ in enumerate(organs):
        node_data = create_cytoscape_node(node_name=organ, node_type='organ', source='miRNATissueAtlas2',
                                          node_data={'id': f'700{idx}', 'display_name': organ})
        network.add_node(organ, **node_data)
    add_edge_from_relationships(network=network, edges=relationships, scores=scores, relationship_type='mo')


def add_organ_system_relationship(network):

    """
    This function will add the organ-system relationships to the network.
    This function will take (for now) the organs as a new node and add the edges to the network.
    The organs are taken from the mirna_tissues table in the database. The relationships are added to the network
    with the add_edge_from_relationships function. The weigths of the edges are set to -1.

    This function will probably be retired.

    :param network: The networkx object
    :return: None

    
    """
    nodes = [x for x, y in network.nodes(data=True) if y['type'] == 'organ']
    relationships, systems = get_tissue_system_edges(nodes)
    for idx, system_n in enumerate(systems):
        node_data = create_cytoscape_node(node_name=system_n, node_type='system', source='miRNATissueAtlas2',
                                          node_data={'id': f'500{idx}', 'display_name': system_n, 'weight': 1})
        network.add_node(system_n, **node_data)
    add_edge_from_relationships(network=network, edges=relationships, relationship_type='os')


def add_other_data(graph, data_df:DataFrame, name:str):
    """

    :param graph:
    :param data_df: Index must be the node (gene name) and the columns the properties.
    :param name: String with the name of that metadata property, examples are cell_type or tissue
    :return: 
    """
    if data_df is None:
        return None
    if name is None:
        return None
    filtered_tissues_df = data_df[data_df.index.isin(graph.nodes)]
    for index, row in filtered_tissues_df.iterrows():
        add_metadata_to_node(graph=graph, node_name=index, tissues=dict(row), name=name)


    pass

def weight_nodes(graph, coeficients: dict):
    """
    This function will add a weigh to each node. This will consist in the sum of the absolute values of the DDS
    times the muscle value

    For now the node's weight is trivial, but this is the  functtion that I will change in 
    https://karenguerrero.atlassian.net/browse/MML-328
    MML-328: Modify the scoring to work dynamically  and not only on one comparison 

    :param graph:
    :return:
    """

    ## Check the data avaialbe in ['data']['metadata']. We
    ## should have dds, pathways, pathways_svd, tf, tissue_expr and cell_type.
    ## Some of them might be empty, in that case, they are ignored. 
    ## There formula is the sum of their absolute values times a coeficient. 
    ## The coeficient is a dictionary with the name of the data as key and the value as the coeficient.
    ## The default value is 1.
    ## The final value is the sum of all the values.
    ## The final value is stored in the node as 'weigh'
    ## A spatial coeficient miR_enhenchment is added if the node is a mirna. It should be 
    ## in the dictionary. If not, the default value is 1.

    if coeficients is None:
        coeficients = {'dds': 1, 'pathways': 1, 'pathways_svd': 1, 'tf': 1, 'tissue_expr': 1, 'cell_type': 1}
    if 'miR_enhancement' not in coeficients:
        coeficients['miR_enhancement'] = 1

    for node, data in graph.nodes(data=True):
        added_coef = []
        weight = 0
        if 'data' in data and 'metadata' in data['data']:

            for key, value in coeficients.items():
                if key in data['data']['metadata']:
                    if key != 'pathways':
                        if key == 'pathways_svd':
                            weight+= value * data['data']['metadata'][key]
                        else:
                            weight += value * sum(abs(x) for x in data['data']['metadata'][key].values())
                        added_coef.append(key)
        if data['type'] == 'mirna':
            weight += coeficients['miR_enhancement']
        graph.nodes[node]['data']['weigh'] = weight






DDS_w =1
yo_w = 1
ym_w = 1
mo_w = 1
miR_w = 1
tissue_w = 1


