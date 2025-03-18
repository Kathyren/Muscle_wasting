import sys
sys.path.append('../../network')

import network.network_processing as nx

pathway_keywords = ["ATP", "MITOCHONDRI", "RESPIRAT", "METABOLI", "OXIDATIVE_PHOSPHORYLATION",
                    "NONALCOHOLIC_FATTY_LIVER", "MUSCLE", "ELECTRON"]
compariosn = ['m_l', 'm_s', 'yo', 'ym', 'mo']

def register_path(graph, node, mir:str, visited_edges=None, path=[]):
    if visited_edges is None:
        visited_edges = []
    node_values = node['data']['influence'][mir]
    node_name = node['data']['name']
    paths= []
    for neighbor in graph.successors(node_name):
        path.append(neighbor)
        edge = graph[node_name][neighbor]
        if edge in visited_edges:
            break
        visited_edges.append(edge)
        edge['weight'] = edge['data']['weight']
        weight = edge['weight']
        if node_name == neighbor:
            effect = node_values[-1]
        else:
            effect = node_values[-1] * weight
        if node_name=="EPAS1" and mir=='hsa-miR-145-5p':
            print('path:', path)
        neighbor_node = graph.nodes[neighbor]
        if 'influence' in neighbor_node['data']:
            if mir in neighbor_node['data']['influence']:
                #print(f"appending_{neighbor}")
                neighbor_node['data']['influence'][mir].append(effect)
            else:
                neighbor_node['data']['influence'][mir] = [effect]
        else:
            neighbor_node['data']['influence'] = {mir: [effect]}
        path = register_path(graph, neighbor_node, mir, visited_edges, path)
        paths.append(path)
    return paths
def visit_all_neighbours(graph, node, mir:str, visited_edges=None):
    if visited_edges is None:
        visited_edges = []

    node_values = node['data']['influence'][mir]
    node_name = node['data']['name']
    for neighbor in graph.successors(node_name):

        edge = graph[node_name][neighbor]
        if edge in visited_edges:
            break
        visited_edges.append(edge)
        if 'data' not in edge:
            print(f"here, {edge}")
        edge['weight'] = edge['data']['weight']
        weight = edge['weight']
        if node_name == neighbor:
            effect = node_values[-1]
        else:
            effect = node_values[-1] * weight
        neighbor_node = graph.nodes[neighbor]
        if 'influence' in neighbor_node['data']:
            if mir in neighbor_node['data']['influence']:
                #print(f"appending_{neighbor}")
                neighbor_node['data']['influence'][mir].append(effect)
            else:
                neighbor_node['data']['influence'][mir] = [effect]
        else:
            neighbor_node['data']['influence'] = {mir: [effect]}
        visit_all_neighbours(graph, neighbor_node, mir, visited_edges)
def start_mir_path(graph, mir):
    node = graph.nodes[mir]
    node['data']['influence']= {mir:[1]}
    visit_all_neighbours(graph, node, mir, visited_edges=[])

def traverse_and_update(graph, start_node):
    # Initialize node values and visited count
    node_values = {node: 0 for node in graph.nodes}
    node_visited = {node: 0 for node in graph.nodes}
    node_marks = {node: 0 for node in graph.nodes}

    def dfs(current_node, current_product):
        for neighbor in graph.successors(current_node):
            print(current_node, neighbor)
            edge = graph[current_node][neighbor]
            if not 'weight' in edge:
                edge['weight'] = edge['data']['weight']

            weight = edge['weight']
            new_product = current_product * weight
            node_visited[neighbor] += 1

            if node_visited[neighbor] > 1:
                if node_values[neighbor] * new_product < 0:
                    node_marks[neighbor] += 1
                node_values[neighbor] += new_product
            else:
                dfs(neighbor, new_product)

    # Start DFS traversal from the start node
    dfs(start_node, 1)

    return node_values, node_marks, node_visited


def join_paths(graph, paths: dict):
    all_nodes = []
    for key, list_paths in paths.items():
        for paths in list_paths:
            all_nodes.extend(paths)

    subgraph = graph.subgraph(all_nodes)
    return subgraph


def evaluate_pathway_influence(influence_mir_data: list, pathway_keywords=None):
    """

    :param pathway_keywords:
    :param influence_data:
    :return:
    """
    if pathway_keywords is None:
        pathway_keywords = ["ATP", "MITOCHONDRI", "RESPIRAT", "METABOLI", "OXIDATIVE_PHOSPHORYLATION",
                            "NONALCOHOLIC_FATTY_LIVER", "MUSCLE", "ELECTRON"]
    pathway_eval_dict = {}
    all_pathways = []
    for pathway in pathway_keywords:
        pathway_repetitions = 0
        for data in influence_mir_data:
            pathway_list = data['pathways']
            pathway_repetitions += sum(1 for s in pathway_list if pathway in s)
            all_pathways.extend(pathway_list)
        pathway_eval_dict[pathway] = pathway_repetitions
    pathway_eval_dict["Different_pathways"] = len(set(all_pathways))
    pathway_eval_dict["Total"] = len(all_pathways)

    return pathway_eval_dict


def evaluate_de_influence(influence_data: list, comparisons):
    if comparisons is None:
        comparisons = compariosn
    comparison_eval_dict = {}
    for comp in comparisons:
        comp_repetitions = 0
        for data in influence_data:
            if data[comp] != 0:
                comp_repetitions += 1
        comparison_eval_dict[comp] = comp_repetitions
    return comparison_eval_dict


def get_influence(graph, paths):
    """
    For each of the paths, visit node by node keeping track of the values on each of the elements
    on metadata.

    :param graph:
    :param paths:
    :return: dictionary with the accumulated influences
    """
    data_list = []
    for path in paths:
        data_dict = evaluate_path_influence(graph=graph, path=path)
        data_list.append(data_dict)
    return data_list


def get_pathways(graph, mirna, n_distance=10, sample_size=100):
    paths = nx.random_walk(graph=graph, node_name=mirna, distance=n_distance, sample_size=sample_size)
    return paths


def evaluate_path_influence(graph, path, comparisons = None):
    """
    This function takes the path, and check the info of the nodes.
    Each node must have the data pathways, yo, ym, mo, ml, ms
    It will check each element and return a dictionary of the final evaluation of the mirna (TBD)

    Goes to the node in the pathm collect the information
    :param graph:
    :param path:
    :return:
    """
    ## Lets take the path and see the nodes

    data_dict = {'weight': 0,
                 'pathways': [],'pathways_svd':0
                 }
    for node in path:
        node_data = graph.nodes[node]
        data_dict = count_influence(node_data=node_data, data_dict=data_dict, comparisons=comparisons)
    return data_dict


def count_influence( node_data, data_dict, comparisons=None):
    # count influence
    if data_dict is None:
        data_dict={'weight':0, 'pathway':[], 'pathways_svd':0}
    if 'data' in node_data:
        if 'weigh' in node_data['data']:
            data_dict['weight'] += node_data['data']['weigh']
        elif 'weight' in node_data['data']:
            data_dict['weight'] += node_data['data']['weigh']
        if 'metadata' in node_data['data']:
            metadata = node_data['data']['metadata']
            if comparisons is None:
                if 'dds' in metadata.keys():
                    comparisons = list(metadata['dds'].keys())
                else:
                    comparisons=[]
            if 'pathways' in metadata:
                pathways = node_data['data']['pathways']
                data_dict['pathways'].extend(pathways)
            if 'pathways_svd' in metadata:
                data_dict['pathways_svd']+= metadata['pathways_svd']

            for comp in comparisons:
                if comp in metadata['dds']:
                    if f"dds_{comp}" in data_dict:
                        data_dict[f"dds_{comp}"] += abs(metadata['dds'][comp])
                    else:
                        data_dict[f"dds_{comp}"] = abs(metadata['dds'][comp])
                if 'tf' in metadata and comp in metadata['tf']:
                    if f"tf_{comp}" in data_dict:
                        data_dict[f"tf_{comp}"] += abs(metadata['tf'][comp])
                    else:
                        data_dict[f"tf_{comp}"] = abs(metadata['tf'][comp])
    return data_dict
