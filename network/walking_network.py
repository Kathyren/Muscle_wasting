import network_processing as nx

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


def evaluate_pathway_influence(influence_data: list):
    """

    :param influence_data:
    :return:
    """
    pathway_eval_dict = {}
    all_pathways = []
    for pathway in pathway_keywords:
        pathway_repetitions = 0
        for data in influence_data:
            pathway_list = data['pathways']
            pathway_repetitions += sum(1 for s in pathway_list if pathway in s)
            all_pathways.extend(pathway_list)
        pathway_eval_dict[pathway] = pathway_repetitions
    pathway_eval_dict["Different_pathways"] = len(set(all_pathways))
    pathway_eval_dict["Total"] = len(all_pathways)

    return pathway_eval_dict


def evaluate_de_influence(influence_data: list):
    comparison_eval_dict = {}
    for comp in compariosn:
        comp_repetitions = 0
        for data in influence_data:
            if data[comp] != 0:
                comp_repetitions += 1
        comparison_eval_dict[comp] = comp_repetitions
    return comparison_eval_dict


def get_influence(graph, paths):
    data_list = []
    for path in paths:
        data_dict = evaluate_path_influence(graph=graph, path=path)
        data_list.append(data_dict)
    return data_list


def get_pathways(graph, mirna, n_distance=10, sample_size=100):
    paths = nx.random_walk(graph=graph, node_name=mirna, distance=n_distance, sample_size=sample_size)
    return paths


def evaluate_path_influence(graph, path):
    """
    This function takes the path, and check the info of the nodes.
    Each node must have the data pathways, yo, ym, mo, ml, ms
    It will check each element and return a dictionary of the final evaluation of the mirna (TBD)
    :param graph:
    :param path:
    :return:
    """
    ## Lets take the path and see the nodes

    data_dict = {'pathways': [],
                 'm_l': 0,
                 'm_s': 0,
                 'yo': 0,
                 'ym': 0,
                 'mo': 0
                 }
    for node in path:
        node_data = graph.nodes[node]
        data_dict = count_influence(graph, node_data, data_dict)
    return data_dict


def count_influence(graph, node_data, data_dict):
    # count influence

    if 'data' in node_data:
        if 'pathways' in node_data['data']:
            pathways = node_data['data']['pathways']
            data_dict['pathways'].extend(pathways)
        for comp in compariosn:
            if comp in node_data['data'] and node_data['data'][comp] != 0:
                data_dict[comp] += 1
    return data_dict
