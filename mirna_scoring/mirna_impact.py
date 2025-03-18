from distutils.core import setup_keywords

import pandas as pd
from scipy.spatial.distance import pdist, squareform
import mirna_scoring.walking_network as wn
import mirna_scoring.jupyter_functions as jf


def count_de_gene(gene_node, comparisons=None, dds_threshold=0):
    """
    for each gene, count how many times it is DE (max is the amount of comparisons)
    :return:
    """
    quantity = 0
    if 'dds' in gene_node['data']['metadata'] :
        dds = gene_node['data']['metadata']['dds']
        for comparison in comparisons:
            comp = dds[comparison]
            # if comp is greater than the threshold, then it will be counted, otherwise it will be 0
            if abs(comp) > dds_threshold:
                quantity += 1
    # quantity is lower 0, max the amount of comparisons
    return quantity


class mirna_network:
    def __init__(self, network):
        self.network = network
        self.miR_nodes = get_mirna_nodes(network)
        self.influence_df = mirnas_influence_on_genes(miR_nodes=self.miR_nodes, network=network) # [[1,1][1,-1]]
        self.influence_sum_df = self.get_impact_data() # [2, 0]
        self.mirnas_paths = None
        #self.mirna_pathway_influence_df = None # get_mirna_pathway_influence_df(network, self.paths)
        self.mirnas_influences = {}
        #self.mirna_de_influence_df = None
        self.pathway_list = self.get_pathway_list()
        self.frequency_pathways = get_frequency_pathways(self.pathway_list)
        #self.top_pathways = None
        self.upRegulated={}
        self.downRegulated={}
        self.available_combinations=self.get_all_available_combinations()
        self.measurements=None
        self.best_mirnas= {}
        self.mirna_distance_matrix=None
        self.mirna_scores={}
        pass
    # ------ Functions to initialized -----
    def get_impact_data(self):
        return get_impact_data(self.influence_df)
    def get_pathway_list(self):
        """
            Get all the pathways of the network
            :param network:
            :return:
            """
        pathway_list = []
        for node in self.network.nodes:
            node = self.network.nodes[node]
            if 'pathways' in node['data']:
                p = node['data']['pathways']
                pathway_list.extend(p)
        return pathway_list
    def get_all_available_combinations(self):
        all_conditions = []
        for node_name in self.network.nodes:
            node = self.network.nodes[node_name]
            if 'data' in node and 'metadata' in node['data']:
                if 'dds' in node['data']['metadata']:
                    all_conditions = node['data']['metadata']['dds'].keys()
        return list(all_conditions)

    # ------ Functions to run -------
    def save_influence_df(self, name):
        self.influence_df.to_csv(name)
    def set_mirnas_paths(self, steps=5,sample_size =10 ):
        """
        This function takes the network and gets for each of the
        mirnas the random walks.
        :param steps:
        :param sample_size:
        :return:
        """
        self.mirnas_paths = get_paths(network=self.network,nodes_start=self.miR_nodes, steps=steps,sample_size=sample_size)
        return self.mirnas_paths
    def set_pathway_influences(self):
        mirna_pathway_influence_df = \
            get_mirna_pathway_influence_df(self.network, self.mirnas_paths)
        self.mirnas_influences['pathway_influences'] = mirna_pathway_influence_df
    def set_de_influences(self):
        mirna_de_influence_df =\
            get_mirna_dde_influence_df(self.network, self.mirnas_paths)
        self.mirnas_influences['DE_influences'] = mirna_de_influence_df
    def set_all_conditions_up_down_regulated(self):
        all_conditions = self.available_combinations
        for condition in all_conditions:
            self.set_up_down_regulated(condition)
        pass
    def set_up_down_regulated(self, condition):
        upRegulated, downRegulated = get_up_down_regulated(self.network, condition)
        self.upRegulated[condition] = upRegulated
        self.downRegulated[condition] = downRegulated
    def get_up_regulated(self, condition):
        if condition in self.upRegulated:
            return self.upRegulated[condition]
        else:
            self.set_up_down_regulated(condition=condition)
        return self.upRegulated[condition]
    def get_down_regulated(self, condition):
        if condition in self.downRegulated:
            return self.downRegulated[condition]
        else:
            self.set_up_down_regulated(condition=condition)
        return self.downRegulated[condition]
    def reset_paths(self, steps=5, sample_size=10):
        self.mirnas_paths = get_paths(network=self.network,
                                      nodes_start=self.miR_nodes,
                                      steps=steps, sample_size=sample_size)
    def get_frequency_pathways(self):
        return self.frequency_pathways
    def get_top_pathways(self):
        return get_most_frequent_pathways(self.frequency_pathways)
    def get_most_impactful_interactions(self, interest_mirna:str):
        order_most_df = get_most_impactful_interactions(self.influence_df, interest_mirna)
        return order_most_df
    def order_length_mirna(self, interest_mirna:str):
        return order_length_mirna(self.influence_df, interest_mirna)
    def order_impact_mirna(self, interest_mirna:str):
        return order_impact_mirna(self.influence_df, interest_mirna)
    def set_mirnas_similarity(self):
        df_numeric = self.influence_sum_df.apply(pd.to_numeric, errors='coerce')
        dist_matrix = pdist(df_numeric, metric='cityblock')
        dist_matrix = squareform(dist_matrix)
        self.mirna_distance_matrix = pd.DataFrame(dist_matrix,
                                                  index=df_numeric.index,
                                                  columns=df_numeric.index)
    def get_mirnas_similarity(self):
        """
        get the mirna distance of influence
        :return:
        """
        if self.mirna_distance_matrix is None:
            self.set_mirnas_similarity()
        return self.mirna_distance_matrix

    def get_mirna_nodes(self):
        return self.miR_nodes
    def get_influences_df(self):
        return self.influence_df
    def get_scores(self):
        return self.mirnas_influences
    def get_best_inhibitor_in_comparison(self, comparison):
        if comparison in self.available_combinations:
           return self.mirnas_influences['inhibitors'][comparison].idxmax()
        else:
            raise Exception(f"No such comparison as {comparison}")
    def get_best_activator_in_comparison(self, comparison):
        if comparison in self.available_combinations:
           return self.mirnas_influences['activators'][comparison].idxmax()
        else:
            raise Exception(f"No such comparison as {comparison}")
    def get_most_pathway_key_word_mirna(self):
        self.mirnas_influences['pathways']['participation'].idxmax()
    # ----- Functions for scoring -----
    def set_measurements(self, comparisons=None, dds_threshold = 0 ):
        if comparisons is None:
            comparisons= self.available_combinations
        self.measurements = self.calculate_measurements(comparisons=comparisons, dds_threshold=dds_threshold)
    def set_up_regulated_dds_score(self):
        """
        Set the mirna influence and calculate the score
        :return: none
        """
        self.mirnas_influences['activators'] = pd.DataFrame(self.calculate_up_regulated_dds_score())
        dds_scores = self.mirnas_influences['activators'].copy()
        # change the name of the columns adding the prefix 'up_'
        dds_scores.columns = ['up_'+col for col in dds_scores.columns]
        # add each combination (column) into a value on mirna_scores, no longer under activators
        for column in dds_scores.columns:
            self.mirna_scores[column] = dds_scores[column]

    def quick_get_all_scores(self, steps=5, sample_size=10, dds_threshold=0, pathway_keywords=[]):
        """

        :param steps:
        :param sample_size:
        :param dds_threshold:
        :param pathway_keywords:
        :return:
        """
        self.set_all_conditions_up_down_regulated()
        self.set_measurements()
        self.set_up_regulated_dds_score()
        self.set_down_regulated_dds_score()
        self.set_quantity_score(dds_threshold=dds_threshold)
        self.set_weight_score()
        self.get_random_walk_pathway_influence(steps=steps,
                                               sample_size=sample_size,
                                               pathway_keywords=pathway_keywords)
        self.set_pathway_score()
        scores = pd.DataFrame(self.mirna_scores)
        return scores



    def set_down_regulated_dds_score(self):
        self.mirnas_influences['inhibitors'] = pd.DataFrame(self.calculate_down_regulated_dds_score())
        dds_scores = self.mirnas_influences['inhibitors'].copy()
        # change the name of the columns adding the prefix 'down_'
        dds_scores.columns = ['down_' + col for col in dds_scores.columns]
        # add each combination (column) into a value on mirna_scores, no longer under inhibitors
        for column in dds_scores.columns:
            self.mirna_scores[column] = dds_scores[column]
    def set_pathway_score(self):
        self.mirnas_influences['pathways']['participation'] = self.mirnas_influences['pathways'].drop(
            columns=["Different_pathways", "Total"]).sum(axis=1)
        pathway_score = self.mirnas_influences['pathways']['participation']
        self.mirna_scores['pathway_counts'] = pathway_score
        self.mirnas_influences['pathway_svd'] = self.calculate_pathway_score()
        self.mirna_scores['pathway_svd'] = self.mirnas_influences['pathway_svd']
    def set_quantity_score(self, dds_threshold=0):
        self.mirnas_influences['de_count'] = self.calculate_quantity_score(dds_threshold)
    def set_weight_score(self):
        weight_score = self.measurements['weight']
        score = weight_score.sum(axis=1)
        self.mirna_scores['gene_weight'] = score


    def calculate_measurements(self, comparisons=None, dds_threshold=0):
        if comparisons is None:
            comparisons = self.available_combinations
        influence_weight = self.influence_sum_df.copy()  # This is the sum of the weights of the nodes it impact
        # multiplied by the times it reaches.
        influence_quantity = self.influence_sum_df.copy()
        influence_pathway = self.influence_sum_df.copy()
        n_mirnas = len(self.influence_sum_df)
        dataframe_dict = {}
        for gene in self.influence_sum_df.columns:
            gene_node = self.network.nodes[gene]
            weight = gene_node['data']['weigh']
            weight_score = weight * self.influence_sum_df[gene]
            influence_weight[gene] = weight_score
            if 'pathways_svd' in gene_node['data']['metadata']:
                pathway_sdv = gene_node['data']['metadata']['pathways_svd']
            else:
                pathway_sdv = 0
            pathway_score = pathway_sdv * self.influence_sum_df[gene]
            influence_pathway[gene] = pathway_score
            quantity = 0
            if 'dds' in gene_node['data']['metadata'] :
                dds = gene_node['data']['metadata']['dds']
                for comparison in comparisons:
                    comp = dds[comparison]
                    # if comp is greater than the threshold, then it will be counted, otherwise it will be 0                        
                    if comparison in dataframe_dict:
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                    else:
                        dataframe_dict[comparison] = {}
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                quantity = count_de_gene(gene_node=gene_node,
                 comparisons=comparisons, dds_threshold=dds_threshold)

            # quantity is lower 0, max the amount of comparisons
            # if abs(self.influence_sum_df[gene]) is >0, then 1, otherwise 0
            affects_gene = self.influence_sum_df[gene].apply(lambda x: 1 if abs(x) > 0 else 0)
            quantity = quantity * affects_gene

            influence_quantity[gene] = abs(quantity)
        for combination, series in dataframe_dict.items():
            series_df = pd.DataFrame(series)
            dataframe_dict[combination] = series_df
        dataframe_dict['weight'] = influence_weight
        dataframe_dict['quantity'] = influence_quantity
        dataframe_dict['pathways'] = influence_pathway

        return dataframe_dict

    def calculate_up_regulated_dds_score(self):
        up_regulated_genes = self.upRegulated
        if self.measurements is None:
            self.measurements = self.calculate_measurements()
        return self.calculate_dds_score(up_regulated_genes)

    def calculate_down_regulated_dds_score(self):
        down_regulated_genes = self.downRegulated
        if self.measurements is None:
            self.measurements = self.calculate_measurements()
        return self.calculate_dds_score(down_regulated_genes)
    def calculate_dds_score(self, up_genes):
        dds_scores= {}
        for comparison in self.available_combinations:
            dds = self.measurements[comparison]
            genes = up_genes[comparison].keys()
            genes = [g for g in genes if g in dds.columns]
            dds = dds[genes]
            dds_score = dds.sum(axis=1)
            dds_scores[comparison] = dds_score
        return dds_scores
    def calculate_pathway_score(self):
        pathway_df = self.measurements['pathways']
        pathway_score = pathway_df.sum(axis=1)
        return pathway_score
    def calculate_quantity_score(self, dds_threshold=0):
        """
        Amount of DE reached (if a gene is DE in two conditions, then it counts for 2)
        :return:
        """
        if 'quantity' not in self.measurements:
            self.calculate_measurements()
        quantity_df = self.measurements['quantity']
        return quantity_df.sum(axis=1)
    def get_random_walk_pathway_influence(self, steps=5, sample_size=10, pathway_keywords=None):
        if self.mirnas_paths is None:
            self.reset_paths(steps=steps, sample_size=sample_size)
        if pathway_keywords is None:
            pathway_keywords = []
        mir_pathway_influence_df = get_mirna_pathway_influence_df(network=self.network,
                                                                  mirPaths=self.mirnas_paths,
                                                                  pathway_keywords=pathway_keywords)
        self.mirnas_influences['pathways']= mir_pathway_influence_df
        return mir_pathway_influence_df













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
        if 'influence' in node['data']:
            influence = node['data']['influence']
            if type(influence) == dict and 'mirna' not in node['type']:
                influences[node['data']['name']] = influence
    df = pd.DataFrame(influences)
    return df

def mirnas_influence_on_genes(miR_nodes, network):
    for mir in miR_nodes:
        try:
            wn.start_mir_path(network, mir)
        except KeyError as ke:
            print(ke)

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
        if 'metadata' in node['data']:
            metadata = node['data']['metadata']
            if 'dds' in metadata and condition in metadata['dds']:
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
        if 'pathways' in node['data']:
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

def something():
    df_impact = jf.get_impact_data(df_down_ym)
    int_influence_df = get_impact_data(influence_df)
    measurements = jf.calculate_measurements(int_influence_df)