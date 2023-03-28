import json

# Load the JSON file into a dictionary
with open('network_from_source_x.json') as f:
    network_data = json.load(f)

# Get the node and edge keys for the source
source = NetworkSource.X
node_keys = source.get_node_keys()
edge_keys = source.get_edge_keys()

# Extract the node data from the JSON file
node_data = []
for node_id, node_attributes in network_data['elements']['nodes'].items():
    node = {key: node_attributes.get(value_key, None) for key, value_key in node_keys.items()}
    node_data.append(node)

# Extract the edge data from the JSON file
edge_data = []
for edge_id, edge_attributes in network_data['elements']['edges'].items():
    edge = {key: edge_attributes.get(value_key, None) for key, value_key in edge_keys.items()}
    edge_data.append(edge)

# Save the extracted data to a file
with open('extracted_data_from_source_x.json', 'w') as f:
    json.dump({'nodes': node_data, 'edges': edge_data}, f)
