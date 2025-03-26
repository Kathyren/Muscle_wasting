
import os

import pandas as pd
import decoupler as dc
## print pwd
import sys
from mirkitten import get_DE_genes_df

from mirkitten.TF import TranFact
from mirkitten.DDS import DDS
import json
import pandas as pd
import random

class new_network:
    def __init__(self, dds_files_path, only_DE=False, threshold=None, pvalue=0.05, interest='stat', species='human', save_name='tf_network.csv'):
        """
        """

        #self.dds_join = pd.read_csv(dds_join)
        self.only_DE = only_DE
        self.de_threshold = threshold
        self.nodes=None
        self.pvalue = pvalue
        self.interest = interest
        dds = DDS(dds_files_path=dds_files_path, pvalue=pvalue, threshold=threshold, interest=interest)
        tf = TranFact(dds_files_path=dds_files_path, pvalue=pvalue, threshold=threshold, interest=interest, species=species)
        self.dds_join = dds.combine_dds_interest_genes(interest=interest)
        self.collectri = tf.colletri
        self.nodes_df = self.collect_genes()
        self.save_name = save_name
        if '.' in self.save_name:
            self.base_name = os.path.basename(self.save_name).split('.')[0]
        else:
            self.base_name = os.path.basename(self.save_name)
        self.source_target_df = None
        pass

    def collect_genes(self):
        """
        This function will check the dataframe with the DDS and
        get all (if only_DE is false) or only the differentially express genes
        :return:
        """
        if self.only_DE:
            de_genes_df = get_DE_genes_df(dds=self.dds_join, pvalue=self.pvalue, threshold=self.de_threshold, interest=self.interest)
            nodes_df = de_genes_df
        else:
            nodes_df = self.dds_join
        self.nodes = nodes_df.index.tolist()
        return nodes_df
    def get_network(self):
        """
        This function will return the network
        """
        selected_genes = self.nodes
        relevant_tf = self.collectri.query('source in @selected_genes or target in @selected_genes')
        #relevant_tf.to_csv(self.save_name)
        self.source_target_df = relevant_tf

    def generate_cytoscape_json(self):
        """
        Generate a Cytoscape JSON file from a DataFrame.
        """
        df = self.source_target_df
        # get the name from save_name (only the name not the directory)
        network_name = self.base_name
        output_file = f"{network_name}.cyjs"
        unique_nodes = set(df['source']).union(set(df['target']))
        node_map = {node: str(i) for i, node in enumerate(unique_nodes, start=1)}  # Assign unique IDs

        nodes = [
            {
                "data": {
                    "id": node_map[node],
                    "shared_name": node,
                    "SUID": int(node_map[node]),
                    "selected": False,
                    "name": node
                },
                "position": {
                    "x": random.uniform(-100, 100),  # Random positions
                    "y": random.uniform(-100, 100)
                },
                "selected": False
            }
            for node in unique_nodes
        ]
        
        edges = [
            {
                "data": {
                    "id": str(i + 1000),  # Ensure unique IDs for edges
                    "source": node_map[row['source']],
                    "target": node_map[row['target']],
                    "shared_name": f"{row['source']} (pm) {row['target']}",
                    "shared_interaction": "pm",
                    "name": f"{row['source']} (pm) {row['target']}",
                    "interaction": "pm",
                    "weight": row['weight'],
                    "SUID": i + 1000,
                    "selected": False
                },
                "selected": False
            }
            for i, row in df.iterrows()
        ]
        
        cytoscape_json = {
            "format_version": "1.0",
            "generated_by": "cytoscape-3.10.1",
            "target_cytoscapejs_version": "~2.1",
            "data": {
                "shared_name": network_name,
                "name": network_name,
                "SUID": 1418,
                "selected": True
            },
            "elements": {
                "nodes": nodes,
                "edges": edges
            }
        }
        
        with open(output_file, "w") as f:
            json.dump(cytoscape_json, f, indent=2)
        
        return cytoscape_json
    
    def save_network(self):
        """
        Save the network as a CSV  and json file
        """
        cytoscape_json = self.generate_cytoscape_json()
        #save the json file
        with open(self.save_name, "w") as f:
            json.dump(cytoscape_json, f, indent=2)
        return self.save_name