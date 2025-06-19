import os
import pandas as pd
import decoupler as dc
## print pwd
import sys
from mirkitten import get_DE_genes_df
import yaml
from mirkitten import get_elbow_point_threshold
from mirkitten.plot_GSEA_ORA import plot_ora_results

import logging
# Set up logging to mirkitten.log
logging.basicConfig(filename='mirkitten.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Pathways:
    def __init__(self, dds_dict, sel_db=None, pvalue=0.05, threshold=None, interest='stat', pathway_pvalue=0.05):
        """
        Initialize the Pathways class with dds_dict and parameters for pathway enrichment analysis.
        :param dds_dict: dict, dictionary with DDS names as keys and DataFrames as values
        :param sel_db: list, list of selected databases for pathway enrichment
        :param pvalue: float, p-value threshold for differential expression
        :param threshold: float, threshold for significance
        :param interest: str, column name in DDS DataFrame to use for significance
        :param pathway_pvalue: float, p-value threshold for pathway enrichment
        """

        logging.info(f'Initializing Pathways with dds_dict: {dds_dict}, sel_db: {sel_db}, pvalue: {pvalue}, threshold: {threshold}, interest: {interest}, pathway_pvalue: {pathway_pvalue}')
        if sel_db is None:
            sel_db = ['go_molecular_function',
                      'go_cellular_component',
                      'go_biological_process',
                      'reactome_pathways',
                      'kegg_pathways', 'hallmark']
        if pathway_pvalue is None:
            pathway_pvalue=0.05
        if "msigdb" in os.listdir('data'):
            msigdb = pd.read_csv('msigdb.csv')
        else:
            msigdb = dc.get_resource('MSigDB')
        self.msigdb =msigdb[msigdb['collection'].isin(sel_db)]
        self.msigdb = self.msigdb[~self.msigdb.duplicated(['geneset', 'genesymbol'])]
        self.dds_dict = self._load_dds_files(dds_files_path=dds_dict)
        self.sel_db = sel_db
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.pathway_pvalue = pathway_pvalue
        self.enriched_ORA_pathway = None
        self.pathways_dict_enriched = {}
    def _load_dds_files(self, dds_files_path) -> dict:
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        :param dds_files_path: str, path to the YAML file with DDS file paths
        :return: dict, dictionary with DDS names as keys and DataFrames as values
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}
    def set_sel_db(self, sel_db):
        self.sel_db = sel_db
    def load_sel_df_from_file(self, path):
        """
        This function will read a text file that has the selected databases
        and take them a list of string that then will be added to self.sel_db
        :param path: str, path to the file with the selected databases
        :return: None
        """
        logging.info(f'Loading selected databases from file: {path}')
        with open(path, 'r') as file:
            self.sel_db = file.readlines()
    
    def select_only_specific_pathways(self,pathways_interest):
        """
        This function will take the pathways of interest
        """
        self.msigdb = self.msigdb[self.msigdb['geneset']==pathways_interest]
        pass
    def restart_msigdb(self)-> None:
        """
        This function will restart the msigdb to the original state.
        It will read the msigdb.csv file from the data folder or will get it from decoupler.
        """
        logging.info('Restarting msigdb to the original state on data/msigdb.csv')
        if "msigdb" in os.listdir('/home/amore/work/data'):
            self.msigdb = pd.read_csv('msigdb.csv')
        else:
            self.msigdb = dc.get_resource('MSigDB')
        self.msigdb =self.msigdb[self.msigdb['collection'].isin(self.sel_db)]
    def show_msigdb_genesets(self):
        return self.msigdb['geneset'].unique()
    def show_msigdb_collections(self):
        return self.msigdb['collection'].unique()
    
    def run_ora(self, degs)-> pd.DataFrame:
        """
        This function will take the DE genes and will run the Pathway enrichment ORA.
        :param degs: list, list of DE genes
        :return: pd.DataFrame, dataframe with the enriched pathways
        """
        logging.info(f'Running ORA with {len(degs)} DE genes')
        self.msigdb = self.msigdb[~self.msigdb.duplicated(['geneset', 'genesymbol'])]
        degs = [item for item in degs if "Unnamed" not in item and "NAN" not in item]
        self.enr_pvals = dc.get_ora_df(
            df=degs,
            net=self.msigdb,
            source='geneset',
            target='genesymbol'
        )
        logging.info(f'Found {len(self.enr_pvals)} enriched pathways')
        return self.enr_pvals

    def get_enriched_pathways(self, df)-> pd.Series:
        """
        This function will takes the DE genes and will run the Pathway enrichment ORA.

        :param df: pd.DataFrame, dataframe with the DE genes
        :return: pd.DataFrame, dataframe with the enriched pathways
        """
        logging.info(f'Getting enriched pathways from DE genes dataframe with shape {df.shape}')
        DE_df = get_DE_genes_df(df, pvalue=self.pvalue, threshold=self.threshold, interest=self.interest)
        degs = DE_df.index
        pathway_df = self.run_ora(degs)
        pathway_df = pathway_df[pathway_df['p-value'] < self.pathway_pvalue]
        if self.threshold is not None:
            pathway_df = pathway_df[pathway_df['Combined score'] > self.threshold]
        else:
            threshold = get_elbow_point_threshold(df = pathway_df, interest='Combined score')
            pathway_df = pathway_df[pathway_df['Combined score'] > threshold]
        enriched_pathways = pathway_df.set_index("Term", inplace=False)  # Set "Term" as index
        enriched_pathways = enriched_pathways['Combined score']
        logging.info(f'Found {len(enriched_pathways)} enriched pathways with p-value < {self.pathway_pvalue} and Combined score > {self.threshold}')
        return enriched_pathways, pathway_df

    def get_enriched_pathways_combined_ORA(self):
        """
        This function will take the dictionary of
        the already DE genes, and will run the Pathway enrichment ORA.
        """
        logging.info('Getting enriched pathways combined ORA')
        pathways_dict = {}
        for comparison, dds in self.dds_dict.items():
            logging.info(f'Processing comparison: {comparison}')
            pathways, all_data = self.get_enriched_pathways(dds)
            self.pathways_dict_enriched[comparison] = all_data
            del all_data
            pathways_dict[comparison] = pathways

        all_scores_pathway_df=None
        for comparison, score in pathways_dict.items():
            if all_scores_pathway_df is None:
                all_scores_pathway_df = score.rename(comparison).to_frame()
            else:
                all_scores_pathway_df = all_scores_pathway_df.join(score.rename(comparison))
            logging.info(f'Added {comparison} to all_scores_pathway_df with shape {all_scores_pathway_df.shape}')
        all_scores_pathway_df.fillna(0, inplace=True)
        return all_scores_pathway_df
    def save_enriched_pathways_combined_ORA(self, path='results/enriched_pathways_combined_ORA.csv'):
        """
        This function will save the enriched pathways combined ORA to a csv file.
        :param path: str, path to save the enriched pathways combined ORA
        """
        logging.info(f'Saving enriched pathways combined ORA to {path}')
        pathway_df = self.get_enriched_pathways_combined_ORA()
        pathway_df.to_csv(path)
    def set_enriched_ORA_pathway(self ):
        """
        This function will set the enriched_ORA_pathway to the enriched pathways combined ORA.
        It will call the get_enriched_pathways_combined_ORA() function.
        """
        logging.info('Setting enriched_ORA_pathway to the enriched pathways combined ORA')
        self.enriched_ORA_pathway = self.get_enriched_pathways_combined_ORA()
        #print (self.enriched_ORA_pathway.head())
        
    def attach_gene_list_to_pathway(self, enrriched_pathway_df=None)-> pd.DataFrame:
        """
        This function will take all the pathways in enriched_ORA_pathway,
          will look at the msigdb and will attach the gene list to the pathway in a new column "genes"
        :param enrriched_pathway_df: pd.DataFrame, dataframe with the enriched pathways
        :return: pd.DataFrame, dataframe with the enriched pathways and the gene list attached

        """
        if enrriched_pathway_df is None:
            if self.enriched_ORA_pathway is None:
                self.set_enriched_ORA_pathway()
            enrriched_pathway_df = self.enriched_ORA_pathway

        # The msigdb has columns genesymbol, collection, geneset.
        # for each pathway in the enriched_ORA_pathway, we will look at the msigdb and filter by the geneset. 
        # It will take all the genes in genesymbol and will attach them to the pathway in a new column "genes" as a list
        msigdb = self.msigdb
        msigdb = msigdb.set_index('geneset')
        gene_list = []

        for pathway in enrriched_pathway_df.index:
            if pathway in msigdb.index:
                genes = msigdb.loc[pathway, "genesymbol"]
                if isinstance(genes, pd.Series):  # If multiple values exist
                    gene_list.append(genes.tolist())
                else:  # Single value case
                    gene_list.append([genes])
                #gene_list.append(genes)
            else: 
                print (pathway)
        enrriched_pathway_df['genes'] = gene_list
        
        return enrriched_pathway_df

    def get_enriched_pathways_with_genes_ORA(self):
        """
        This function will return the enriched pathways with the genes
        """
        if self.enriched_ORA_pathway is None:
            self.set_enriched_ORA_pathway()
        return self.attach_gene_list_to_pathway()
    
    def save_enrriched_pathways_with_genes_ORA(self, path):
        """
        This function will save the enriched pathways with genes to a csv file.
        :param path: str, path to save the enriched pathways with genes
        """
        logging.info(f'Saving enriched pathways with genes to {path}')
        df = self.get_enriched_pathways_with_genes_ORA()
        df.to_csv(path)

    def plot_oras(self, combination, filename, title="", top_n=10, figsize=(12, 6), scale_odds_ratio=.5,
                     fontsize_title=12, fontsize_subtitle=12, fontsize_text=10):
        """
        This function will plot the ORA results for a given combination of pathways.
        :param combination: str, the name of the combination of pathways to plot
        :param filename: str, the name of the file to save the plot
        :param title: str, the title of the plot
        :param top_n: int, the number of top pathways to plot
        :param figsize: tuple, the size of the figure
        :param scale_odds_ratio: float, the scale of the odds ratio
        :param fontsize_title: int, the font size of the title
        :param fontsize_subtitle: int, the font size of the subtitle
        :param fontsize_text: int, the font size of the text
        """

        if self.pathways_dict_enriched == {}:
            raise ValueError('You need to run the get_enriched_pathways_combined_ORA() function first')
        pathway_df = self.pathways_dict_enriched[combination]
        if len(pathway_df) == 0:
            logging.warning(f'No pathways found for combination {combination}. Skipping plot.')
            return
        fig, ax = plot_ora_results(pathway_df, top_n=top_n, figsize=figsize, scale_odds_ratio=scale_odds_ratio,
                                  fontsize_title=fontsize_title, fontsize_subtitle=fontsize_subtitle,
                                  fontsize_text=fontsize_text, title=title)
        fig.savefig(filename)
        logging.info(f'Saved ORA plot to {filename}')
    def save_individual_pathway_results(self, save_individual, n_top=10, plots=False):
        """
        This function will save the individual pathway results to a file.
        :param save_individual: str, the path to save the individual results
        :param n_top: int, the number of top pathways to save
        :param plots: bool, whether to save the plots or not
        """
        logging.info(f'Saving individual pathway results to {save_individual} with top {n_top} pathways')
        if os.path.isfile(save_individual):
            save_individual = os.path.dirname(save_individual)
        if save_individual[-1] != '/':
            save_individual += '/'
        for comparison in self.pathways_dict_enriched.keys():
            pathway_df = self.pathways_dict_enriched[comparison]
            logging.info(f'Processing comparison: {comparison}, type of pathway_df: {type(pathway_df)}')
            logging.info(f'Processing comparison: {comparison} with shape {pathway_df.shape}')
            pathway_df = pathway_df.sort_values(by='Combined score', ascending=False)
            pathway_df.to_csv(os.path.join(save_individual, f'pathway_de_{comparison}.csv'))
            if plots:
                self.plot_oras(combination=comparison,
                               filename=os.path.join(save_individual, f'pathway_de_{comparison}.svg'),
                               title=f'Pathway enrichment for {comparison}', top_n=n_top)
            logging.info(f'Saved individual pathway results for {comparison} to {save_individual}pathway_de_{comparison}.csv')


