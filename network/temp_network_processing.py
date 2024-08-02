import random
import networkx as nx
import numpy as np
import pandas as pd
from database_analysis import sql_operations as sql
from cytoscape import format_cytoscape_json, create_cytoscape_node, create_cytoscape_edge, protein_name
import pickle


class MiRNAGeneNetwork:
    def __init__(self):
        self.graph = None

    def create_graph(self, mirnas, genes, relationship):
        """
        This function will receive a list of genes and mirnas and create the network of them
        :param relationship:
        :param mirnas: List of string with the mirna names
        :param genes: List of string with the gene names
        :param relationsip: List of tuple of string with the node names
        """
        self.graph = nx.Graph()
        self.graph.add_nodes_from(mirnas, size=10, type="mirna")
        self.graph.add_nodes_from(genes, type="gene")
        self.graph.add_edges_from(relationship)
        return self.graph

    def create_graph_from_dictionaries(self, nodes, edges):
        """
        This function will receive a list of genes and mirnas and create the network of them
        :param edges: List of dictionaries with the metadata of each relationship
        :param nodes: List of dictionary with the node metadata
        """

        self.graph = nx.DiGraph()

        for element in nodes:
            node = element['data'][protein_name]
            if node not in self.graph._node:
                self.graph.add_node(node, **element)
        for edge in edges:
            try:
                if edge['target'] not in self.graph._adj[edge['source']]:
                    self.graph.add_edge(edge['source'], edge['target'], **edge)
            except Exception as e:
                print(e)
        return self.graph

    def draw_graph(self):
        color_map = []
        for node, data in self.graph.nodes(data=True):
            if data['type'] == 'mirna':
                color_map.append(0.25)  # blue color
            elif data['type'] == 'gene':
                color_map.append(0.7)  # yellow color
        nx.draw(self.graph, node_color=color_map, with_labels=False)

    def filter_by_degree(self):
        graph = self.graph.copy()
        centrality = nx.degree(self.graph)
        c = list(centrality)
        _, values = zip(*c)
        average = np.mean(values)
        standard_dev = np.std(values)
        x = np.quantile(values, 0.75)
        delete_nodes = [node[0] for node in centrality if node[1] < x]
        graph.remove_nodes_from(delete_nodes)
        self.save_graph(graph, "high_degree.pkl")
        return graph

    def filter_by_threshold(self, pathway_file, pathway_feature, threshold_feature, threshold_value):
        pathway_df = pd.read_csv(pathway_file, index_col=0)
        pathways = pathway_df.loc[pathway_df[threshold_feature] > threshold_value, pathway_feature]
        genes = set()
        for path in pathways:
            genes.update(path.split(','))
        return genes

    def remove_nodes_low_centrality(self, cutoff=0.75):
        cc_t = nx.closeness_centrality(self.graph)
        cc = list(cc_t.values())
        x = np.quantile(cc, cutoff)
        delete_nodes = [node for node, value in cc_t.items() if value < x]
        delete_edges = [(u, v) for node in delete_nodes for u, v in self.graph.edges(node)]
        self.graph.remove_edges_from(delete_edges)
        self.graph.remove_nodes_from(delete_nodes)
        return self.graph

    def evaluate_edge_quality(self):
        required_keys = {"id", "interaction", "name", "selected", "shared_interaction", "shared_name", "source",
                         "target", "weight"}
        nodes_with_ids = {data['data']['id']: node for node, data in self.graph.nodes(data=True) if
                          'data' in data and isinstance(data['data'], dict) and 'id' in data['data']}
        reverse_dict = {v: k for k, v in nodes_with_ids.items()}

        edges_to_remove = []

        for u, v, data in self.graph.edges(data=True):
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

        self.graph.remove_edges_from(edges_to_remove)

    def evaluate_node_quality(self):
        required_keys = {"id", "SUID", "name", "selected", "shared_name"}

        nodes_to_remove = []

        for node, data in self.graph.nodes(data=True):
            if 'data' in data:
                if isinstance(data['data'], dict):
                    missing_keys = required_keys - data['data'].keys()
                    for key in missing_keys:
                        data['data'][key] = None  # or any default value you prefer
                else:
                    nodes_to_remove.append(node)
            else:
                nodes_to_remove.append(node)

        self.graph.remove_nodes_from(nodes_to_remove)

    def remove_node_and_edges(self, node):
        delete_nodes = [node]
        delete_edges = []
        delete_edges.extend(self.graph.in_edges(node))
        delete_edges.extend(self.graph.out_edges(node))
        self.graph.remove_edges_from(delete_edges)
        return self.graph

    def mark_TF_nodes_from_file(self, tf_file):
        tf_dic = pd.read_csv(tf_file, index_col=0).to_dict('index')

        for node, data in self.graph.nodes(data=True):
            if 'data' in data and 'id' in data['data'] and data['data']['name'] in tf_dic.keys():
                self.graph.nodes[node]['data']['node_type'] = 'TF'
                self.graph.nodes[node]['data']['presence'] = tf_dic[data['data']['name']]['enrichment']

    def mark_miR_nodes(self):
        for node, data in self.graph.nodes(data=True):
            if str('data' in data and 'id' in data['data'] and data['data']['name']).startswith('hsa-miR'):
                self.graph.nodes[node]['data']['node_type'] = 'miR'
                self.graph.nodes[node]['data']['presence'] = None

    def get_n_mirs(self):
        mir_nodes = 0
        for node, data in self.graph.nodes(data=True):
            if 'data' in data and 'node_type' in data['data'] and data['data']['node_type'] == 'miR':
                mir_nodes += 1

        return mir_nodes

    def get_mirs(self):
        mir_nodes = []
        for node, data in self.graph.nodes(data=True):
            if 'data' in data and 'node_type' in data['data'] and data['data']['node_type'] == 'miR':
                mir_nodes.append(node)

        return mir_nodes


    def add_DDS_data(self, dds_df):
        for index, row in dds_df.iterrows():
            self.add_dds_to_node(node_name=index, dds=dict(row))

    def add_dds_to_node(self, node_name, dds: dict):
        for node, data in self.graph.nodes(data=True):
            if 'data' in data and 'id' in data['data'] and data['data']['name'] == node_name:
                for dd, value in dds.items():
                    self.graph.nodes[node]['data'][dd] = value

    def add_tissue_to_node(self, node_name, tissues: dict):
        for node, data in self.graph.nodes(data=True):
            if 'data' in data and 'id' in data['data'] and data['data']['name'] == node_name:
                for index, value in tissues.items():
                    tissues[index] = float(value)
                self.graph.nodes[node]['data']['tissue_expr'] = str(tissues)

    def add_tissue_to_nodes(self, tissues_df):
        filtered_tissues_df = tissues_df[tissues_df.index.isin(self.graph.nodes)]
        for index, row in filtered_tissues_df.iterrows():
            self.add_tissue_to_node(node_name=index, tissues=dict(row))

    def save_graph(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self.graph, f)

    def load_graph(self, filename):
        with open(filename, 'rb') as f:
            self.graph = pickle.load(f)
        return self.graph

    def get_interest_genes_and_neighbors(self, n_neighbors: int, distance: int, interest_genes: list):
        """
        This function will check a list of genes of interest and only keep those genes, and their nodes at distance n_neighbors
        :param n_neighbors:
        :param distance:
        :param interest_genes:
        :return:
        """
        for gene in interest_genes:
            pass
        pass

    def select_next_valid_node(self, i, path_length, neighbors) -> str:
        next_node = random.choice(neighbors)
        if i < (path_length - 1):
            n_neighbors = len(list(self.graph.successors(next_node)))
            while len(neighbors) > 1 > n_neighbors:
                data = self.graph.nodes[next_node]
                if 'data' in data and 'pathways' in data['data']:
                    break
                neighbors.remove(next_node)
                next_node = random.choice(neighbors)
                n_neighbors = len(list(self.graph.successors(next_node)))
        return next_node

    def get_random_path(self, start_node, path_length):
        path = [start_node]
        current_node = start_node

        for i in range(path_length - 1):
            try:
                neighbors = list(
                    self.graph.successors(current_node))  # Get the successors (out-neighbors) of the current node

            except KeyError as ke:
                return []
            if not neighbors:
                break  # No more neighbors to explore, end the path here

            next_node = self.select_next_valid_node(i, path_length, neighbors)

            path.append(next_node)
            current_node = next_node

        return path

    def random_walk(self, node_name, distance, sample_size) -> list:
        """
        This function takes the graph, looks for the corresponding node, and does the random path.
        :param node_name:
        :param distance:
        :param sample_size:
        :return:
        """
        if node_name is None:
            return []
        paths = []
        for i in range(sample_size):
            path = self.get_random_path(node_name, distance)
            paths.append(path)
        return paths

    def remove_nodes_low_centrality_pageRank(self, weight=None, cutoff=0.75):
        """
        This function will calculate a centrality

        :param cutoff: What is going to be the quartile of nodes that we are going to keep
        :param weight: The weight parameter for the PageRank algorithm
        :return: The graph without the un-relevant nodes
        """
        cc_t = nx.pagerank(self.graph, weight=weight)
        cc = list(cc_t.values())
        x = np.quantile(cc, cutoff)
        delete_nodes = []
        delete_edges = []
        for node in cc_t.items():
            if node[1] < x:
                delete_nodes.append(node[0])
                delete_edges.extend(self.graph.edges(node[0]))
        self.graph.remove_edges_from(delete_edges)
        self.graph.remove_nodes_from(delete_nodes)
        return self.graph

    def convert_to_json(self):
        """
        This function will take a networkx object to json
        :return: json object
        """
        json_data = nx.cytoscape_data(self.graph)
        json_data = nx.node_link_data(self.graph)
        json_data['elements'] = {'nodes': json_data.pop('nodes'), 'edges': json_data.pop('links')}
        return json_data

    def get_nodes_names(self):
        """
        Return the name of the nodes in a list of string
        :return: list of node names
        """
        return list(self.graph.nodes)

    def get_mirna_mrna_relationships(self, genes):
        """
        This function will take the gene names and get the list
        of microrna and their relationships

        :param genes: list of string; List of gene names
        :return: The list of tuples for the relationships and the list of string with the name of the mirnas and the scores
        """
        genes_str = '"' + '","'.join(genes) + '"'
        query = f'Select Distinct mrna, mirna_mature from gene_mirna where mirna_mature like "hsa-%" and mrna in ({genes_str});'
        relationships = sql.get_query(query=query)
        genes = list(relationships['mrna'])
        mirnas = list(relationships['mirna_mature'])
        relationship = list(zip(mirnas, genes))
        scores = [1] * len(genes)
        return relationship, mirnas, scores

    def get_mirna_tissue_edges(self, mirnas):
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

    def get_tissue_system_edges(self, organs):
        """
        This function will collect the organ-system edges
        :return: tuple of the relationships and the systems
        """
        organ_list = '"' + '","'.join(organs) + '"'
        query = f"select distinct mirna_tissues.system, organ from mirna_tissues where organ in ({organ_list});"
        result = sql.get_query(query=query)
        organ = list(result['organ'])
        systems = list(result['system'])
        relationship = list(zip(organ, systems))
        return relationship, systems

    def create_my_graph(self, query):
        relationships = sql.get_query(query=query)
        genes = list(relationships['mrna'])
        mirnas = list(relationships['mirna_mature'])
        relationship = list(zip(genes, mirnas))
        self.graph = self.create_graph(mirnas=mirnas, genes=genes, relationship=relationship)
        return self.graph

    def add_relationships(self, relationships):
        """
        This function will add to an existent network the relationships
        :param relationships: A list of tuples with the relationships
        :return: None
        """
        self.graph.add_edges_from(relationships)

    def set_positions(self):
        """
        This function will rearrange the nodes and update the position as specified by cytoscape
        :return:
        """
        self.evaluate_node_quality()
        pos = nx.spectral_layout(self.graph, scale=100, center=[50, 50])

        for node in self.graph.nodes:
            try:
                self.graph.nodes[node]['position']['x'] = pos[node][0]
                self.graph.nodes[node]['position']['y'] = pos[node][1]
            except Exception as e:
                print(f"Error with node {node}: {str(e)}")
                self.remove_node_and_edges(node)

        self.evaluate_edge_quality()


    def add_mirna_relationships(self):
        genes = self.get_nodes_names()
        relationships, mirnas, scores = self.get_mirna_mrna_relationships(genes)
        for idx, mirna in enumerate(mirnas):
            node_data = create_cytoscape_node(node_name=mirna, node_type='mirna', source='mirbase',
                                              node_data={'id': f'900{idx}', 'display_name': mirna, 'shared_name': mirna,
                                                         'name': mirna})
            self.graph.add_node(mirna, **node_data)
        self.add_edge_from_relationships(edges=relationships, scores=scores)

    def add_edge_from_relationships(self, edges, scores=None, relationship_type="pm"):
        if scores is None:
            scores = [1] * len(edges)
        for idx, edge in enumerate(edges):
            source, target = edge
            try:

                if 'data' not in self.graph.nodes[source]:
                    self.graph.nodes[source]['data'] = {'id': source}
                if 'data' not in self.graph.nodes[target]:
                    self.graph.nodes[target]['data'] = {'id': target}
                source_data = self.graph.nodes[source]['data']
                target_data = self.graph.nodes[target]['data']
                data2fill = {"id": f'800{idx}', "source": source_data['id'], "target": target_data['id'],
                             "shared_interaction": relationship_type,
                             "shared_name": f"{source} ({relationship_type}) {target}",
                             "name": f"{source} ({relationship_type}) {target}", "interaction": relationship_type,
                             "selected": False, "weight": -1}
                edge_data = create_cytoscape_edge(source=source, target=target, node_data=data2fill, weight=scores[idx])
                self.graph.add_edge(source, target, **edge_data)
            except Exception as e:
                print(f"Relationship {source}-{target} couldn't be added due to {str(e)}")
        self.add_relationships(edges)

    def add_tissue_relationship(self):
        nodes = self.get_nodes_names()
        relationships, organs, scores = self.get_mirna_tissue_edges(nodes)
        for idx, organ in enumerate(organs):
            node_data = create_cytoscape_node(node_name=organ, node_type='organ', source='miRNATissueAtlas2',
                                              node_data={'id': f'700{idx}', 'display_name': organ})
            self.graph.add_node(organ, **node_data)
        self.add_edge_from_relationships(edges=relationships, scores=scores, relationship_type='mo')

    def add_organ_system_relationship(self):
        nodes = [x for x, y in self.graph.nodes(data=True) if y.get('type') == 'organ']
        relationships, systems = self.get_tissue_system_edges(nodes)
        for idx, system_n in enumerate(systems):
            node_data = create_cytoscape_node(node_name=system_n, node_type='system', source='miRNATissueAtlas2',
                                              node_data={'id': f'500{idx}', 'display_name': system_n, 'weight': 1})
            self.graph.add_node(system_n, **node_data)
        self.add_edge_from_relationships(edges=relationships, relationship_type='os')

    def weight_nodes(self, tissues=None):
        """
        This function will add a weight to each node. This will consist in the sum of the absolute values of the DDS
        times the muscle value
        :param tissues:
        :return:
        """
        if tissues is None:
            tissues = ['skeletal muscle']
        for node, data in self.graph.nodes(data=True):
            if 'data' not in data:
                continue
            total_w = 0
            if 'yo' in data['data']:
                total_w += abs(data['data']['yo'])
            if 'mo' in data['data']:
                total_w += abs(data['data']['mo'])
            if 'ym' in data['data']:
                total_w += abs(data['data']['ym'])
            if 'tissue_expr' in data['data']:
                for tissue in tissues:
                    total_w += (eval(data['data']['tissue_expr'])[tissue]) / len(tissues)
            self.graph.nodes[node]['data']['weigh'] = total_w + 0.01

        for node, data in self.graph.nodes(data=True):
            if 'data' not in data:
                continue
            if 'node_type' in data['data'] and data['data']['node_type'] == 'miR':
                total_w = 0
                all_neighbors = nx.all_neighbors(self.graph, node)
                n = 0
                for neighbor in all_neighbors:
                    neighbor_node = self.graph.nodes[neighbor]
                    if 'data' in neighbor_node and 'weigh' in neighbor_node['data']:
                        neighbor_w = neighbor_node['data']['weigh']
                        n += 1
                    else:
                        neighbor_w = 0
                    total_w += neighbor_w
                self.graph.nodes[node]['data']['weigh'] = total_w - 0.01 * (n - 1)


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


def extract_genes_from_pathways(pathway_file, threshold_feature='Combined score',
                                id_feature='Features', pathway_feature='Term', threshold_value=50):
    """
    This function takes a fiile that cointains
    id term set size overlap ratio p_value fdr p-vale odds ratio combine score and features


    :param threshold_feature:
    :param id_feature:
    :param pathway_feature:
    :param pathway_file:
    :param threshold_value:
    :return:
    """
    pathway_df = pd.read_csv(pathway_file)
    pathway_df = pathway_df[pathway_df[threshold_feature] > threshold_value]
    feature_dict = {}

    # Iterate through each row in the dataframe
    for index, row in pathway_df.iterrows():
        # Split the Features column by ';' to get individual features
        features = row[id_feature].split(';')

        # Iterate through each feature
        for feature in features:
            if feature:  # Ignore empty strings
                if feature not in feature_dict:
                    feature_dict[feature] = []
                # Append the Term (pathway) to the list of pathways for this feature
                feature_dict[feature].append(row[pathway_feature])
    return feature_dict


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


get_random_relationships = "Select Distinct mrna, mirna_mature from binding limit 10;"
get_muscle_relationships = 'Select Distinct mrna, mirna_mature from binding where  ' \
                           'mirna_mature like "hsa-%" and mrna in (select gene from gene_bank) ;'
