import pandas as pd
from scipy.spatial.distance import pdist, squareform
import mirna_scoring.walking_network as wn



def get_impact_data(df):
    """

    :param df: Dataframe of the influences with the [-1,1] etc
    :return:
    """
    def convert_to_int_list(lst):
        if isinstance(lst, list):
            return [int(x) for x in lst]
        else:
            return []

    color_data = pd.DataFrame(index=df.index, columns=df.columns)

    for col in df.columns:
        for idx in df.index:
            int_list = convert_to_int_list(df.at[idx, col])
            color_data.at[idx, col] = sum(int_list)

    return color_data

def order_length_mirna(df, interest_mirna:str):
    df_sorted = df.copy()
    df_sorted['length'] = df_sorted[interest_mirna].apply(lambda x: len(x) if isinstance(x, list) else 0)
    df_sorted = df_sorted.sort_values(by='length', ascending=False)
    df_sorted = df_sorted.drop(columns=['length'])
    return df_sorted


def order_impact_mirna(df, interest_mirna:str):
    df_sorted = df.copy()
    df_sorted['impact'] = df_sorted[interest_mirna].apply(lambda x: sum(x) if isinstance(x, list) else 0)
    df_sorted = df_sorted.sort_values(by='impact', ascending=False)
    df_sorted = df_sorted.drop(columns=['impact'])
    return df_sorted

def get_most_impactful_interactions(df, interest_mirna):
    order_df = order_impact_mirna(df, interest_mirna)
    highest = order_df[:100]
    lowest = order_df[-100:]
    order_most_df = pd.concat([highest, lowest])
    return order_most_df

def get_mirnas_similarity(df, interest_mirna): #todo: This function is not done
    """
    This funtion will return the matrix of distance of the effect of pairs of microRNAs
    :return:
    """

    impact_most = get_impact_data(order_most_df)

    df_numeric = impact_most.apply(pd.to_numeric, errors='coerce')

    # Calculate the pairwise distances between columns
    dist_matrix = pdist(df_numeric.T, metric='cityblock')
    dist_matrix_square = squareform(dist_matrix)

    # Create a DataFrame for the distance matrix
    dist_df = pd.DataFrame(dist_matrix_square, index=df_numeric.columns, columns=df_numeric.columns)

    return dist_df

def sum_lengths(row):
    return sum(len(lst) for lst in row)

# Add the new column 'X'
#todo: the way I fins mirnas is not correct. use node['type'] instead
def get_mirna_nodes(network):
    miR_nodes = [node for node in network.nodes if node.startswith('hsa-miR')]
    return miR_nodes

def get_influences_df(network):
    """
    From the network as input, get the property 'influence' and get a dataframe
    with the genes in the columns and the mirnas in the rows.
    :param network:
    :return:
    """
    influences = {}
    for node_name in network.nodes:
        node = network.nodes[node_name]
        if 'data' in node and 'influence' in node['data']:
            influence = node['data']['influence']
            if type(influence) == dict and 'mirna' not in node['type']:
                influences[node['data']['name']] = influence
        else:
            pass
            #print(f"Node {node_name} does not have data or influence, it has {node}")
    df = pd.DataFrame(influences)
    return df

def mirnas_influence_on_genes(miR_nodes:list, network):
    for mir in miR_nodes:
       # try:
            wn.start_mir_path(network, mir)
        #except KeyError as ke:
        #    print(f"Error: mirna {mir} does not have {ke}")

    influences = get_influences_df(network)
    influences = influences[influences.index.isin(miR_nodes)]
    influences = influences.fillna("").apply(list)

    return influences

def get_up_down_regulated(network, condition):
    """

    :param network:
    :param condition:
    :return:
    """
    regulation_YO = {}
    up_regulation_YO = {}
    down_regulation_YO = {}
    for node_name in network.nodes:
        node = network.nodes[node_name]
        if 'data' in node and 'metadata' in node['data']:
            metadata = node['data']['metadata']
            if 'dds_original' in metadata and condition in metadata['dds_original']:
                regulation = metadata['dds_original'][condition]
                regulation_YO[node_name] = regulation
                if regulation > 0:
                    up_regulation_YO[node_name] = regulation
                elif regulation < 0:
                    down_regulation_YO[node_name] = regulation
            elif 'dds' in metadata and condition in metadata['dds']:
                regulation = metadata['dds'][condition]
                regulation_YO[node_name] = regulation
                if regulation > 0:
                    up_regulation_YO[node_name] = regulation
                elif regulation < 0:
                    down_regulation_YO[node_name] = regulation
    return up_regulation_YO, down_regulation_YO

def get_up_down_regulated_df(df, network, condition):
    """
    This funtion will take the already influence dataframe with genes in the columns
    and mirnas in the rows
    :param df:
    :param network:
    :param condition:
    :return:
    """
    up_regulation_YO, down_regulation_YO = get_up_down_regulated(network, condition)
    df_down_yo = df[df.index.isin(down_regulation_YO.keys())]
    df_up_yo = df[df.index.isin(up_regulation_YO.keys())]
    return df_up_yo, df_down_yo
def get_paths(network, nodes_start:list, steps:int=5, sample_size=10):
    """
    This function will run for each node in node_start (idealy mirnas) the random walk in the network and
    get what are the nodes it went thorough
    :param network: The network to be analized
    :param nodes_start: The list of nodes (mirnas) that are going to be the path start for the random walk
    :param steps:
    :param sample_size:
    :return: Dictionary with the node_starts and their respective paths
    """
    steps = steps + 1
    mirPaths = {}
    for mir in nodes_start:
        p = wn.get_pathways(graph=network, mirna=mir, n_distance=steps, sample_size=sample_size)
        unique_set = set(tuple(lst) for lst in p)

        # Convert back to a list of lists
        unique_list = [list(tpl) for tpl in unique_set]
        mirPaths[mir] = unique_list

    return mirPaths

def get_mirna_influence(mirPaths:dict, network):
    """
    For each mirna in mirPaths, get what is their accumulated influence
    :param mirPaths:
    :param network:
    :return: dictionary with the mirs and the accumulated influence
    """
    mirInfluence = {}
    for mir, path in mirPaths.items():
        influence = wn.get_influence(network, path)
        mirInfluence[mir] = influence
    return mirInfluence
def get_mirna_pathway_influence_df(network, mirPaths, pathway_keywords=None):
    """
    This function will take the mirPath calculated previously
    and just sum how many times it landed on a gene with a specific pathway
    It will return a dataframe with the pathways in the rows and how many times it landed on a specific
    pathway.
    :param network:
    :param mirPaths:
    :return:
    """
    mirInfluence = get_mirna_influence(mirPaths=mirPaths, network=network)
    mir_pathway_influence = {}
    for mir, influence_data in mirInfluence.items():
        pi = wn.evaluate_pathway_influence(influence_mir_data=influence_data, pathway_keywords=pathway_keywords)
        mir_pathway_influence[mir] = pi
    mir_pathway_influence_df = pd.DataFrame(mir_pathway_influence).T
    return mir_pathway_influence_df

def get_mirna_dde_influence_df(network, mirPaths ):
    mirInfluence = get_mirna_influence(mirPaths=mirPaths, network=network)
    mir_de_influence = {}
    for mir, influence_data in mirInfluence.items():
        pi = wn.evaluate_de_influence(influence_data)
        mir_de_influence[mir] = pi
    mir_de_influence_df = pd.DataFrame(mir_de_influence).T
    return mir_de_influence_df
def get_mirnas_with_de_count(mir_de_influence_df, min_val = 0 ):
    mirs_all_infliences = mir_de_influence_df[
        (mir_de_influence_df['ym'] > min_val) |
        (mir_de_influence_df['mo'] > min_val) |
        (mir_de_influence_df['yo'] > min_val)]
    return mirs_all_infliences

def get_pathway_list(network):
    """
    Get all the pathways of the network
    :param network:
    :return:
    """
    pathway_list = []
    for node in network.nodes:
        if 'data' in node and 'pathways' in node['data']:
            p = node['data']['pathways']
            pathway_list.extend(p)
    return pathway_list

def get_frequency_pathways(pathway_list):
    frequency_pathways = dict((i, pathway_list.count(i)) for i in pathway_list)
    frequency_pathways = dict(sorted(frequency_pathways.items(), key=lambda item: item[1]))
    return frequency_pathways

def get_most_frequent_pathways(frequency_pathways, n=20):
    top = {}
    for key, value in frequency_pathways.items():
        if value > n:
            top[key] = value
    return top

