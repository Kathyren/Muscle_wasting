import networkx as nx
import numpy as np

import network.node_evaluation
from database_analysis import sql_operations as sql
from cytoscape import format_cytoscape_json, create_cytoscape_node, create_cytoscape_edge, protein_name


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
            if edge['target'] not in G._adj[edge['source']] :
                G.add_edge(edge['source'], edge['target'], **edge)
        except Exception as e:
            print(e)

    # is_bipartite(G)

    return G


def draw_graph(G):
    color_map = []

    for node, data in G.nodes(data=True):
        if data['type'] == 'mirna':
            color_map.append(0.25)  # blue color
        elif data['type'] == 'gene':
            color_map.append(0.7)  # yellow color
    nx.draw(G, node_color=color_map, with_labels=False)
    # nx.draw(G)


def is_bipartite(G):
    return True
    for edge, data in G.edges(data=True):
        count = 0
        for node in edge:
            if '-miR-' in node:
                count = count + 1
        assert count == 1, 'There is a mirna_mirna interaction or a gene- gene interaction'
    return True


def filter_by_degree(G):
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
    delete_edges=[]
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

def remove_nodes_low_centrality_pageRank(graph, cutoff=0.75):
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
    cc_t = nx.pagerank(graph)
    cc = list(cc_t.values())
    x = np.quantile(cc, cutoff)
    delete_nodes = []
    delete_edges=[]
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
    relationship = list(zip(genes, mirnas))
    scores = [1]*len(genes) # list(relationships['probability'])

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
    genes = get_nodes_names(network)
    relationships, mirnas, scores = get_mirna_mrna_relationships(genes)
    for idx, mirna in enumerate(mirnas):
        node_data = create_cytoscape_node(node_name=mirna, node_type='mirna', source='mirbase',
                                          node_data={'id': f'900{idx}', 'display_name': mirna, 'shared_name':mirna, 'name':mirna})
        network.add_node(mirna, **node_data)

    add_edge_from_relationships(network=network, edges=relationships, scores=scores)
    # add_relationships(network, relationships=relationships)


def add_edge_from_relationships(network, edges, scores=None, relationship_type="pm"):
    # edges, scores = get_mirna_tissue_edges(mirnas=mirnas)
    if scores is None:
        scores= [1]*len(edges)
    for idx, edge in enumerate(edges):
        try:
            source = edge[0]
            target = edge[1]
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
                         "weight": scores[idx]}
            edge_data = create_cytoscape_edge(source=source, target=target, node_data=data2fill, weight=scores[idx])
            network.add_edge(source, target, **edge_data)
        except Exception as e:
            print(f"Relationship {source}-{target} couldn't be added due to {str(e)}")
    add_relationships(network, edges)


def add_tissue_relationship(network):
    nodes = get_nodes_names(network)
    relationships, organs, scores = get_mirna_tissue_edges(nodes)
    for idx, organ in enumerate(organs):
        node_data = create_cytoscape_node(node_name=organ, node_type='organ', source='miRNATissueAtlas2',
                                          node_data={'id': f'700{idx}', 'display_name': organ})
        network.add_node(organ, **node_data)
    add_edge_from_relationships(network=network, edges=relationships, scores=scores, relationship_type='mo')


def add_organ_system_relationship(network):
    nodes = [x for x, y in network.nodes(data=True) if y['type'] == 'organ']
    relationships, systems = get_tissue_system_edges(nodes)
    for idx, system_n in enumerate(systems):
        node_data = create_cytoscape_node(node_name=system_n, node_type='system', source='miRNATissueAtlas2',
                                          node_data={'id': f'500{idx}', 'display_name': system_n, 'weight': 1})
        network.add_node(system_n, **node_data)
    add_edge_from_relationships(network=network, edges=relationships, relationship_type='os')


if __name__ == '__main__':
    # graph = load_graph()
    c_nodes, c_edges, c_relationships = format_cytoscape_json()
    graph = create_graph_from_dictionaries(nodes=c_nodes, relationship=c_relationships, edges=c_edges)
    # graph = create_my_graph(get_random_relationships)
    save_graph(graph, "small_graph.pkl")
    anc = nx.all_pairs_node_connectivity(graph)

    print(nx.info(graph))

    # draw_graph(graph)
