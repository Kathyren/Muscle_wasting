import os
import pandas as pd
import decoupler as dc
import yaml
import json
import logging
import mirkitten.plot_TF as plot_TF
# Set up logging to mirkitten.log
logging.basicConfig(filename='mirkitten.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class TranFact():
    def __init__(self, dds_files_path, pvalue=0.05, threshold=None, interest='stat', species='human'):
        logging.info(f'Initializing TranFact with dds_files_path: {dds_files_path}, pvalue: {pvalue}, threshold: {threshold}, interest: {interest}, species: {species}')
        dds_dict = self._load_dds_files(dds_files_path=dds_files_path)
        self.dds_dict = dds_dict
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.colletri = self.set_colletri(species=species)
        self.mats = {}
        self.tf_dfs = {}
        self.plotter = None
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}
    
    def set_colletri(self, species='human'):
        if species=='human' and "collectri.csv" in os.listdir('data'):
            logging.info('Loading genes for human, reading collectri from data/collectri.csv')
            collectri = pd.read_csv('data/collectri.csv', index_col=0)
        else:
            logging.info(f'Loading genes for {species}, getting collectri from decoupler')
            collectri = dc.get_collectri(organism=species, split_complexes=True)
            logging.info(f'Collection ready. Collectri for {species} has {len(collectri)} genes')
        self.collectri = collectri
        return collectri

    def get_mat(self, dds, interest='stat', comparison='comparison'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        logging.info(f'Getting matrix for comparison: {comparison} with interest: {interest}')
        if comparison in self.mats:
            mat = self.mats[comparison]
        if f'mat_{comparison}_{interest}.csv' in os.listdir('results'):
            mat = pd.read_csv(f'results/mat_{comparison}.csv', index_col=0)
            mat = mat.T
        else:
            mat = dds[[interest]].T.rename(index={interest: comparison})
            
            try:
                mat.T.to_csv(f'results/mat_{comparison}_{interest}.csv')
            except OSError:
                print(f"coudn't save at results/mat_{comparison}_{interest}.csv. Currently at {os.system('pwd')}")
        self.mats[comparison]= mat
        return mat
    def get_tf_df(self, mat, comparison='comparison'):
        logging.info(f'Getting transcription factor dataframe for comparison: {comparison}')
        if comparison in self.tf_dfs and self.tf_dfs[comparison] is not None:
            tf_df = self.tf_dfs[comparison]
        elif f'tf_acts_{comparison}.csv' in os.listdir('results'):
            tf_df = pd.read_csv(f'results/tf_acts_{comparison}_{self.interest}.csv', index_col=0)
            self.tf_dfs[comparison] = tf_df
        else:
            tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=self.collectri, verbose=True)
            tf_df = pd.DataFrame(tf_acts.T)
            tf_df['pvals']=tf_pvals.T
            try:
                tf_df.to_csv(f'results/tf_acts_{comparison}_{self.interest}.csv')
            except OSError:
                print(f"coudn't save at rresults/tf_acts_{comparison}_{self.interest}.csv. Currently at {os.system('pwd')}")
            self.tf_dfs[comparison] = tf_df
        logging.info(f'Transcription factor dataframe for comparison: {comparison} has {len(tf_df)} genes')
        return tf_df

    def get_tf_df_of_comparison(self, dds, interest='stat', comparison='comparison'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        logging.info(f'Getting transcription factor dataframe for comparison: {comparison} with interest: {interest}')
        mat = self.get_mat(dds=dds, interest=interest, comparison=comparison)
        tf_df = self.get_tf_df(mat=mat, comparison=comparison)
        return tf_df
    def get_tf_df_of_combined_comparisons(self):
        logging.info('Getting transcription factor dataframe for all comparisons combined')
        all_tf_dic = {}
        for comparison, dds in self.dds_dict.items():
            tf_df = self.get_tf_df_of_comparison(dds=dds, comparison=comparison, interest=self.interest)
            tf_df = tf_df[tf_df['pvals']<self.pvalue]
            all_tf_dic[comparison] = tf_df[comparison]
            #print( tf_df.head())
        all_scores_tf_df=None
        for comparison, score in all_tf_dic.items():
            if all_scores_tf_df is None:
                all_scores_tf_df = score.rename(comparison).to_frame()
            else:
                all_scores_tf_df = all_scores_tf_df.join(score.rename(comparison))
        all_scores_tf_df.fillna(0, inplace=True)
        logging.info(f'Combined transcription factor dataframe has {len(all_scores_tf_df)} genes across all comparisons')
        return all_scores_tf_df
    
    def save_tf_df_of_combined_comparisons(self, path='results/tf_df_combined.csv'):
        logging.info(f'Saving transcription factor dataframe for all comparisons combined to {path}')
        tf_df = self.get_tf_df_of_combined_comparisons()
        tf_df.to_csv(path)
        
    def get_enrtiched_tf(self, comparison='comparison', pvalue=0.1):
        """
        This function will return the transcription factors that are enriched in the comparison
        """
        logging.info(f'Getting enriched transcription factors for comparison: {comparison} with pvalue: {pvalue}')
        if self.tf_dfs[comparison] is None:
            raise ValueError('You need to run the get_tf_df() function first')
        tf_df = self.tf_dfs[comparison]
        tf_df = tf_df[tf_df['pvals']<pvalue]
        logging.info(f'Found {len(tf_df)} enriched transcription factors for comparison: {comparison}')
        return tf_df
    
    def get_highes_lowest_tf(self, comparison='comparison', top=5) -> list:
        """
        This function will return the top up and down regulated transcription factors
        """
        logging.info(f'Getting highest and lowest transcription factors for comparison: {comparison}')
        if top is None:
            top = 5
        if self.tf_dfs[comparison] is None:
            raise ValueError('You need to run the get_tf_df() function first')
        tf_df = self.get_enrtiched_tf(comparison=comparison)
        tf_acts=pd.DataFrame(tf_df[comparison])
        tf_acts=tf_acts.T
        values = tf_acts.iloc[0]
        down_reg = values.sort_values(ascending=True)[:top].index.to_list()
        up_reg = values.sort_values(ascending=False)[:top].index.to_list()
        logging.info(f'Found {len(down_reg)} down-regulated and {len(up_reg)} up-regulated transcription factors for comparison: {comparison}')
        up_down_reg = down_reg.copy()
        up_down_reg.extend(up_reg)
        return up_down_reg

    def save_individual_tf_results(self, save_individual, n_top=5, plots=False):
        # check if it is a path or a file
        if os.path.isfile(save_individual):
            save_individual = os.path.dirname(save_individual)
        if save_individual[-1] != '/':
            save_individual += '/'
        tf_up_down_top = {}
        for comparison, tf in self.tf_dfs.items():
            tf_de_df = self.get_enrtiched_tf(comparison=comparison, pvalue=self.pvalue)
            tf_up_down = self.get_highes_lowest_tf(comparison=comparison, top=n_top)
            tf_de_df.to_csv(os.path.join(save_individual, f'tf_de_{comparison}.csv'))
            tf_up_down_top[comparison] = tf_up_down
            if plots:
                self.plot_tf_bar(comparison=comparison, filename=os.path.join(save_individual, f'tf_bar_{comparison}.svg'), top=n_top)
            # Save the dictionary as a json file
        file_name = os.path.join(save_individual, f'tf_up_down_top_{n_top}.json')
        with open(file_name, 'w') as f:
            json.dump(tf_up_down_top, f)
        
    
    def plot_tf_bar(self, comparison, filename, top=5, title_fontsize=14, figsize=(5, 5), vertical=True, dpi=700):
        """
        This function will plot the transcription factors bar plot
        """
        logging.info(f'Plotting transcription factors bar plot for comparison: {comparison}')
        plotter = plot_TF.plot_TF(comparison=comparison, tf_df=self.tf_dfs[comparison])
        plotter.plot_personalize_bars(title=f"Top enriched TF for {comparison}", top=top,
                            sig = 'padj', title_fontsize=title_fontsize, 
                            figsize=figsize, vertical=vertical, dpi=dpi, filename=filename)
        logging.info(f'Saved transcription factors bar plot to {filename}')