import pandas as pd
import numpy as np
import os
## print pwd
os.system('pwd')
import sys
sys.path.append('/home/amore/work')
from mirkitten import get_DE_dds_dict
import yaml
from mirkitten import get_elbow_point_threshold, get_DE_genes_df
import logging
import itertools
import json
# Set up logging to mirkitten.log
logging.basicConfig(filename='mirkitten.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class DDS:
    """
    Class to prepare the DDS input of mirKat.

    It will take the individual DDSs of the different comparisons and will combine them into a single DataFrame.
    
    
    """
    def __init__(self, dds_files_path='../data/dds_files.yml', pvalue=0.05, threshold=0, interest='stat'):
        """
        Initializes the DDS class by loading the DDS files specified in a YAML file.
        """
        self.dds_dict = self._load_dds_files(dds_files_path)
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.dds_df = self.combine_dds_interest_genes(interest=interest)
    
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}

    def combine_dds_interest_genes(self, interest='stat'):
        """
        Combines statistical values of interest from DDS comparisons into a single DataFrame.
        """    
        interest_dict = {}
        for comparison, dds_sig in self.dds_dict.items():
            stat = dds_sig[interest]
            interest_dict[comparison] = stat
        all_stat_df=None
        for comparison, stat in interest_dict.items():
            #comparison = comparison[:-5]
            if all_stat_df is None:
                all_stat_df = stat.rename(comparison).to_frame()
            else:
                all_stat_df = all_stat_df.join(stat.rename(comparison))
        all_stat_df.fillna(0, inplace=True)
        return all_stat_df
    
    def get_dds_df(self):
        """
        Returns the combined DDS DataFrame.
        """
        return self.dds_df
    def save_dds_df(self, path='results/dds_df.csv'):
        """
        Saves the combined DDS DataFrame to a CSV file.
        """
        self.dds_df.to_csv(path)
    def get_DE_genes(self, dds, pvalue=0.05):
        """
        Returns a DataFrame of differentially expressed genes based on the combined DDS DataFrame.
        """
        
        if self.threshold == 0:
            threshold = get_elbow_point_threshold(df=dds, interest=self.interest)
        else:
            threshold = self.threshold
        logging.info(f'Using threshold: {threshold} for interest: {self.interest}')
        logging.info(f"The columns of the dds are: {dds.columns}")
        de_genes_df = get_DE_genes_df(dds, pvalue=pvalue, threshold=threshold, interest=self.interest)
        return de_genes_df
    
    def save_all_DE_genes_df(self, path='results/', pvalue=0.05):
        """
        Saves the DataFrame of differentially expressed genes to a CSV file.
        """
        degs = {}
        for comparison, dds in self.dds_dict.items():
            de_genes_df = self.get_DE_genes(dds, pvalue=pvalue)
            de_genes_df.to_csv(os.path.join(path, f'{comparison}_DE_genes.csv'))
            degs[comparison] = de_genes_df.index.tolist()
        degs = self.get_intersections_dict(degs)
        filename = os.path.join(path, 'intersections_DE_genes.json')
        with open(filename, 'w') as f:
            json.dump(degs, f, indent=2)

    
    def get_intersections_dict(self, gene_sets_dict, min_genes=1):
        """
        Given a dictionary where keys are condition names and values are lists or sets of genes,
        return a dictionary of intersections,
        including only those with more than `min_genes` genes.
        
        Parameters:
            gene_sets_dict (dict): Keys are group names, values are lists/sets of genes.
            min_genes (int): Minimum number of genes in intersection to include.
        
        Returns:
            dict: Keys are concatenated group names, values are lists of intersecting genes.
        """
        intersections_dict = {}
        set_names = list(gene_sets_dict.keys())

        for r in range(2, len(set_names) + 1):
            for combination in itertools.combinations(set_names, r):
                intersection = set(gene_sets_dict[combination[0]])
                for set_name in combination[1:]:
                    intersection &= set(gene_sets_dict[set_name])
                if len(intersection) > min_genes:
                    key = "_".join(combination)
                    intersections_dict[key] = list(intersection)
        logging.info(f'Found {len(intersections_dict)} intersections with more than {min_genes} genes.')
        return intersections_dict
    