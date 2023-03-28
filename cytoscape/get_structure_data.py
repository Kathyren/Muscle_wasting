import json

# Load the JSON file into a dictionary
def get_source_data(network_file_name, source, main_name):
    with open(network_file_name) as f:
        network_data = json.load(f)
    if "elements" in network_data and "nodes" in network_data['elements'] and "edges" in network_data['elements']:
        # Get the node and edge keys for the source
        node_data = list(network_data['elements']['nodes'][0]["data"].keys())
        edge_data = list(network_data['elements']['edges'][0]["data"].keys())
        add_source_data_to_file(source=source, node_data=node_data, edge_data=edge_data, main_name=main_name)
    else:
        print(f"The file presented does not have the correct format.")


    # Extract the node data from the JSON file
def get_node_data(network_data, node_keys):
    node_data = []
    for node_id, node_attributes in network_data['elements']['nodes']:
        node = {key: node_attributes.get(value_key, None) for key, value_key in node_keys}
        node_data.append(node)
    return node_data

    # Extract the edge data from the JSON file
def get_edge_data(network_data, edge_keys):
    edge_data = []
    for edge_id, edge_attributes in network_data['elements']['edges'].items():
        edge = {key: edge_attributes.get(value_key, None) for key, value_key in edge_keys}
        edge_data.append(edge)
    return edge_data

def add_source_data_to_file(source, node_data, edge_data, main_name, file_name="Data/Apps_details.json"):
    # Save the extracted data to a file
    with open(file_name, 'a') as f:
        data = {'source':source, "node_keys":node_data, "edge_keys":edge_data, "main_name": main_name }
        f.write(f'"{source}:"'+json.dumps(data)+"\n")



get_source_data(network_file_name="../network/Networks_CYJS/mo_GeneMania_Base.cyjs", source="gene_mania", main_name="gene_name")