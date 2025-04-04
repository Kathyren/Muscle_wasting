
from mirna_scoring.score import *
from mirna_scoring.mirna_impact_helper import *
import mirna_scoring.jupyter_functions as jf
import mirna_scoring.mirna_influence_plots as mi_plot
from itertools import chain
from fractions import Fraction
class mirna_network:
    def __init__(self, network):
        self.network = network
        self.miR_nodes = get_mirna_nodes(network)
        self.influence_df = mirnas_influence_on_genes(miR_nodes=self.miR_nodes, network=network) # [[1,1][1,-1]]
        self.influence_df_paths=None
        self.influence_sum_df = self.get_impact_data() # [2, 0]
        self.influence_sum_df_paths = None
        self.mirnas_paths = None
        #self.mirna_pathway_influence_df = None # get_mirna_pathway_influence_df(network, self.paths)
        self.mirnas_influences = {}
        #self.mirna_de_influence_df = None
        self.pathway_list = self.get_pathway_list()
        self.frequency_pathways = get_frequency_pathways(self.pathway_list)
        #self.top_pathways = None
        self.upRegulated={}
        self.downRegulated={}
        self.available_combinations=self.get_all_available_combinations()
        self.measurements=None
        self.best_mirnas= {}
        self.mirna_distance_matrix=None
        self.mirna_scores={}
        pass
    # ------ Functions to initialized -----
    def get_impact_data(self):
        return get_impact_data(self.influence_df)
    def get_pathway_list(self):
        """
            Get all the pathways of the network
            :param network:
            :return:
            """
        pathway_list = []
        for node in self.network.nodes:
            node = self.network.nodes[node]
            if 'data' in node and 'pathways' in node['data']:
                p = node['data']['pathways']
                pathway_list.extend(p)
        return pathway_list
    def get_all_available_combinations(self):
        all_conditions = []
        for node_name in self.network.nodes:
            node = self.network.nodes[node_name]
            if 'data' in node and 'metadata' in node['data']:
                if 'dds' in node['data']['metadata']:
                    all_conditions = node['data']['metadata']['dds'].keys()
        return list(all_conditions)

    # ------ Functions to run -------
    def save_influence_df(self, name):
        self.influence_df.to_csv(name)
    def set_mirnas_paths(self, steps=5,sample_size =10, seed=42):
        """
        This function takes the network and gets for each of the
        mirnas the random walks.
        :param seed:
        :param steps:
        :param sample_size:
        :return:
        """
        self.mirnas_paths = get_paths(network=self.network,nodes_start=self.miR_nodes,
                                      steps=steps,sample_size=sample_size, seed=seed)
        return self.mirnas_paths
    def set_pathway_influences(self):
        mirna_pathway_influence_df = \
            get_mirna_pathway_influence_df(self.network, self.mirnas_paths)
        self.mirnas_influences['pathway_influences'] = mirna_pathway_influence_df

    def set_genes_influences(self):
        influences = {}
        for mirna, paths in self.mirnas_paths.items():
            nodes = list(chain.from_iterable(paths))
            sub_network = self.network.subgraph(nodes)
            influences[mirna] = mirnas_influence_on_genes(miR_nodes=self.mirnas_paths.keys(), network=sub_network)
        return influences
    def set_de_influences(self):
        mirna_de_influence_df =\
            get_mirna_dde_influence_df(self.network, self.mirnas_paths)
        self.mirnas_influences['DE_influences'] = mirna_de_influence_df
    def set_all_conditions_up_down_regulated(self):
        all_conditions = self.available_combinations
        for condition in all_conditions:
            self.set_up_down_regulated(condition)
        pass
    def set_up_down_regulated(self, condition):
        upRegulated, downRegulated = get_up_down_regulated(self.network, condition)
        self.upRegulated[condition] = upRegulated
        self.downRegulated[condition] = downRegulated
    def get_up_regulated(self, condition):
        if condition in self.upRegulated:
            return self.upRegulated[condition]
        else:
            self.set_up_down_regulated(condition=condition)
        return self.upRegulated[condition]
    def get_down_regulated(self, condition):
        if condition in self.downRegulated:
            return self.downRegulated[condition]
        else:
            self.set_up_down_regulated(condition=condition)
        return self.downRegulated[condition]
    def reset_paths(self, steps=5, sample_size=10, seed=42):
        self.mirnas_paths = get_paths(network=self.network,
                                      nodes_start=self.miR_nodes,
                                      steps=steps, sample_size=sample_size, seed=seed)
    def get_frequency_pathways(self):
        return self.frequency_pathways
    def get_top_pathways(self):
        return get_most_frequent_pathways(self.frequency_pathways)
    def get_most_impactful_interactions(self, interest_mirna:str):
        order_most_df = get_most_impactful_interactions(self.influence_df, interest_mirna)
        return order_most_df
    def order_length_mirna(self, interest_mirna:str):
        return order_length_mirna(self.influence_df, interest_mirna)
    def order_impact_mirna(self, interest_mirna:str):
        return order_impact_mirna(self.influence_df, interest_mirna)
    def set_mirnas_similarity(self):
        df_numeric = self.influence_sum_df.apply(pd.to_numeric, errors='coerce')
        dist_matrix = pdist(df_numeric, metric='cityblock')
        dist_matrix = squareform(dist_matrix)
        self.mirna_distance_matrix = pd.DataFrame(dist_matrix,
                                                  index=df_numeric.index,
                                                  columns=df_numeric.index)
    def get_mirnas_similarity(self):
        """
        get the mirna distance of influence
        :return:
        """
        if self.mirna_distance_matrix is None:
            self.set_mirnas_similarity()
        return self.mirna_distance_matrix

    def get_mirna_nodes(self):
        return self.miR_nodes
    def get_influences_df(self):
        return self.influence_df
    def get_scores(self):
        return self.mirnas_influences
    def get_best_inhibitor_in_comparison(self, comparison):
        if comparison in self.available_combinations:
           return self.mirnas_influences['inhibitors'][comparison].idxmax()
        else:
            raise Exception(f"No such comparison as {comparison}")
    def get_best_activator_in_comparison(self, comparison):
        if comparison in self.available_combinations:
           return self.mirnas_influences['activators'][comparison].idxmax()
        else:
            raise Exception(f"No such comparison as {comparison}")
    def get_most_pathway_key_word_mirna(self):
        self.mirnas_influences['pathways']['participation'].idxmax()
    # ----- Functions for scoring -----
    def set_measurements(self, comparisons=None, dds_threshold = 0 ):
        if comparisons is None:
            comparisons= self.available_combinations
        self.measurements = self.calculate_measurements(comparisons=comparisons, dds_threshold=dds_threshold)
    def set_up_regulated_dds_score(self):
        """
        Set the mirna influence and calculate the score
        :return: none
        """
        self.mirnas_influences['activators'] = pd.DataFrame(self.calculate_up_regulated_dds_score())
        dds_scores = self.mirnas_influences['activators'].copy()
        # change the name of the columns adding the prefix 'up_'
        dds_scores.columns = ['up_'+col for col in dds_scores.columns]
        # add each combination (column) into a value on mirna_scores, no longer under activators
        for column in dds_scores.columns:
            self.mirna_scores[column] = dds_scores[column]

    def quick_get_all_scores(self, steps=5, sample_size=10, dds_threshold=0, pathway_keywords=[]):
        """

        :param steps:
        :param sample_size:
        :param dds_threshold:
        :param pathway_keywords:
        :return:
        """
        self.set_all_conditions_up_down_regulated()
        self.set_measurements(dds_threshold=dds_threshold)
        self.set_up_regulated_dds_score()
        self.set_down_regulated_dds_score()
        self.set_quantity_score(dds_threshold=dds_threshold)
        self.set_weight_score()
        self.get_random_walk_pathway_influence(steps=steps,
                                               sample_size=sample_size,
                                               pathway_keywords=pathway_keywords)
        self.set_pathway_score()
        scores = pd.DataFrame(self.mirna_scores)
        return scores



    def set_down_regulated_dds_score(self):
        self.mirnas_influences['inhibitors'] = pd.DataFrame(self.calculate_down_regulated_dds_score())
        dds_scores = self.mirnas_influences['inhibitors'].copy()
        # change the name of the columns adding the prefix 'down_'
        dds_scores.columns = ['down_' + col for col in dds_scores.columns]
        # add each combination (column) into a value on mirna_scores, no longer under inhibitors
        for column in dds_scores.columns:
            self.mirna_scores[column] = dds_scores[column]
    def set_pathway_score(self):
        self.mirnas_influences['pathways']['participation'] = self.mirnas_influences['pathways'].drop(
            columns=["Different_pathways", "Total"]).sum(axis=1)
        pathway_score = self.mirnas_influences['pathways']['participation']
        self.mirna_scores['pathway_counts'] = pathway_score
        self.mirnas_influences['pathway_svd'] = self.calculate_pathway_score()
        self.mirna_scores['pathway_svd'] = self.mirnas_influences['pathway_svd']
    def set_quantity_score(self, dds_threshold=0):
        self.mirnas_influences['de_count'] = self.calculate_quantity_score(dds_threshold)
        self.mirna_scores['de_count'] = self.mirnas_influences['de_count']
    def set_weight_score(self):
        weight_score = self.measurements['weight']
        score = weight_score.sum(axis=1)
        self.mirna_scores['gene_weight'] = score


    def calculate_measurements(self, comparisons=None, dds_threshold=0):
        if comparisons is None:
            comparisons = self.available_combinations
        influence_weight = self.influence_sum_df.copy()  # This is the sum of the weights of the nodes it impact
        # multiplied by the times it reaches.
        influence_quantity = self.influence_sum_df.copy()
        influence_pathway = self.influence_sum_df.copy()
        n_mirnas = len(self.influence_sum_df)
        dataframe_dict = {}
        for gene in self.influence_sum_df.columns:
            gene_node = self.network.nodes[gene]
            weight = gene_node['data']['weigh']
            weight_score = weight * self.influence_sum_df[gene]
            influence_weight[gene] = weight_score
            if 'pathways_svd' in gene_node['data']['metadata']:
                pathway_sdv = gene_node['data']['metadata']['pathways_svd']
            else:
                pathway_sdv = 0
            pathway_score = pathway_sdv * abs(self.influence_sum_df[gene])
            influence_pathway[gene] = pathway_score
            quantity = 0
            if 'dds_original' in gene_node['data']['metadata'] :
                dds = gene_node['data']['metadata']['dds_original']
                for comparison in comparisons:
                    comp = dds[comparison]
                    # if comp is greater than the threshold, then it will be counted, otherwise it will be 0                        
                    if comparison in dataframe_dict:
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                    else:
                        dataframe_dict[comparison] = {}
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                quantity = count_de_gene(gene_node=gene_node,
                 comparisons=comparisons, dds_threshold=dds_threshold)
            elif 'dds' in gene_node['data']['metadata'] :
                dds = gene_node['data']['metadata']['dds']
                for comparison in comparisons:
                    comp = dds[comparison]
                    # if comp is greater than the threshold, then it will be counted, otherwise it will be 0
                    if comparison in dataframe_dict:
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                    else:
                        dataframe_dict[comparison] = {}
                        dataframe_dict[comparison][gene] = comp * self.influence_sum_df[gene]
                quantity = count_de_gene(gene_node=gene_node,
                 comparisons=comparisons, dds_threshold=dds_threshold)


            # quantity is lower 0, max the amount of comparisons
            # if abs(self.influence_sum_df[gene]) is >0, then 1, otherwise 0
            affects_gene = self.influence_sum_df[gene].apply(lambda x: 1 if abs(x) > 0 else 0)
            quantity = quantity * affects_gene

            influence_quantity[gene] = abs(quantity)
        for combination, series in dataframe_dict.items():
            series_df = pd.DataFrame(series)
            dataframe_dict[combination] = series_df
        dataframe_dict['weight'] = influence_weight
        dataframe_dict['quantity'] = influence_quantity
        dataframe_dict['pathways'] = influence_pathway

        return dataframe_dict

    def calculate_up_regulated_dds_score(self):
        up_regulated_genes = self.upRegulated
        if self.measurements is None:
            self.measurements = self.calculate_measurements()
        return self.calculate_dds_score(up_regulated_genes)

    def calculate_down_regulated_dds_score(self):
        down_regulated_genes = self.downRegulated
        if self.measurements is None:
            self.measurements = self.calculate_measurements()
        return self.calculate_dds_score(down_regulated_genes)
    def calculate_dds_score(self, up_genes):
        dds_scores= {}
        for comparison in self.available_combinations:
            dds = self.measurements[comparison]
            genes = up_genes[comparison].keys()
            genes = [g for g in genes if g in dds.columns]
            dds = dds[genes]
            dds_score = dds.sum(axis=1)
            dds_scores[comparison] = dds_score
        return dds_scores
    def calculate_pathway_score(self):
        pathway_df = self.measurements['pathways']
        pathway_score = pathway_df.sum(axis=1)
        return pathway_score
    def calculate_quantity_score(self, dds_threshold=0):
        """
        Amount of DE reached (if a gene is DE in two conditions, then it counts for 2)
        :return:
        """
        if 'quantity' not in self.measurements:
            self.calculate_measurements(dds_threshold=dds_threshold)
        quantity_df = self.measurements['quantity']
        return quantity_df.sum(axis=1)
    def get_random_walk_pathway_influence(self, steps=5, sample_size=10, pathway_keywords=None):
        if self.mirnas_paths is None:
            self.reset_paths(steps=steps, sample_size=sample_size)
        if pathway_keywords is None:
            pathway_keywords = []
        mir_pathway_influence_df = get_mirna_pathway_influence_df(network=self.network,
                                                                  mirPaths=self.mirnas_paths,
                                                                  pathway_keywords=pathway_keywords)
        self.mirnas_influences['pathways']= mir_pathway_influence_df
        return mir_pathway_influence_df
    def get_minra_influence_as_fraction(self, mirna, gene):
        # Get the float value
        float_values = self.influence_df.loc[mirna][gene]
        fractions = []
        for float_value in float_values:
            fraction = Fraction(float_value).limit_denominator(100)
            fraction_str = f"{fraction.numerator}/{fraction.denominator}"
            fractions.append(fraction_str)


class mirna_evaluation:

    def __init__(self,mirna_network:mirna_network,
                session_name="Session_1",
                base_directory= 'mirna_scoring/results/' ):
        self.minas_cluster = None
        self.mirna_network = mirna_network
        self.mirna_scores = None
        if base_directory[-1]!='/':
            base_directory+='/'
        self.directory = base_directory+session_name
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            print("Directory created successfully!")
        else:
            print("Directory already exists!")
        self.selected_mirnas = []
        self.dist_df = self.mirna_network.get_mirnas_similarity()
        self.pathway_df = None
        self.clustered_mirnas = None
        self.set_pathway_database()

    def get_all_mirnas(self):
        return self.mirna_network.miR_nodes
    def score(self, steps=10, sample_size=10, dds_threshold=2, pathway_keywords=None):
        if pathway_keywords is None:
            pathway_keywords = []
        self.scores= self.mirna_network.quick_get_all_scores(steps=steps, sample_size=sample_size,
                                                             dds_threshold=dds_threshold,
                                                             pathway_keywords=pathway_keywords)
        return self.scores
    def select_mirnas(self):
        if self.scores is None:
            raise Exception("There is need to score the mirnas before select them")
        self.selected_mirnas=[]
        top_rows_extreme, row_scores_extreme = select_extreme_rows(self.scores, x=5, method='iqr')
        top_rows_ranked, row_scores_top = select_top_ranked_rows(self.scores, x=5)
        top_rows_ranked_normalized, row_scores_normalized = select_top_normalized_rows(self.scores, x=5)

        self.selected_mirnas.extend(top_rows_extreme.index)
        self.selected_mirnas.extend(top_rows_ranked.index)
        self.selected_mirnas.extend(top_rows_ranked_normalized.index)

        for comparison in self.mirna_network.get_all_available_combinations():
            mir = self.mirna_network.get_best_inhibitor_in_comparison(comparison)
            self.selected_mirnas.append(mir)
            mir = self.mirna_network.get_best_activator_in_comparison(comparison)
            self.selected_mirnas.append(mir)
        return self.selected_mirnas

    def get_mirna_influence(self, mirna):
        """
        Get the genes and their influence
        :return:
        """
        paths = self.mirna_network.mirnas_paths[mirna]
        sub_network = []
        for path in paths:
            sub_network.extend(path)
        impact = self.mirna_network.influence_sum_df.loc[mirna]
        return sub_network, impact

    def cluster_mirnas(self, n_clusters=None):
        unique_mirnas = list(set(self.selected_mirnas))
        if n_clusters is None:
            n_clusters = len(unique_mirnas) // 3
        mirna_clusters = jf.cluster_mirnas(dist_matrix_square=self.dist_df, n_clusters=n_clusters)
        minas_cluster = mirna_clusters.sort_values(by=["Cluster"])['Cluster']
        clustered_mirnas = {}
        for mirna in unique_mirnas:
            if mirna in minas_cluster.index:
                cluster = minas_cluster[mirna]
                if cluster not in clustered_mirnas:
                    cluster = int(cluster)
                    clustered_mirnas[cluster] = []
                clustered_mirnas[cluster].append(mirna)
        self.clustered_mirnas = clustered_mirnas
        self.minas_cluster = minas_cluster
        return clustered_mirnas
    def set_pathway_database(self, sel_db=None):
        if sel_db is None:
            sel_db = ['go_molecular_function',
                      'go_cellular_component',
                      'go_biological_process',
                      'reactome_pathways',
                      'kegg_pathways', 'hallmark']
        if "msigdb.csv" in os.listdir('mirkitten/data'):
            self.msigdb = pd.read_csv('mirkitten/data/msigdb.csv', index_col=0)
        else:
            self.msigdb = dc.get_resource('MSigDB')
        self.msigdb.index = self.msigdb['genesymbol']

        msigdb = self.msigdb[self.msigdb['collection'].isin(sel_db)]
        self.msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
    def create_report_single_mirnas(self):
        genes = []
        for mirna in self.selected_mirnas:
            g = self.create_report_single_mirna(mirna)
            genes.extend(g)
        return genes
    def create_report_single_mirna(self, mirna):
        genes = []
        sub_network, impact = self.get_mirna_influence(mirna)

        sub_network = list(set(sub_network))
        # mi_plot.draw_network(G=my_network.network, node_list = sub_network, name= mirna, save_path='mirna_scoring/sub_plots')
        html_file = mi_plot.generate_html_network_report(network=self.mirna_network.network, selected_nodes=sub_network,
                                                         name=mirna, save_path=self.directory)
        # fig, ax = get_plot_enriched(sub_network, f"Enriched pathways of targets of {mirna}")
        fig, ax = get_mirna_target_enriched(selected_genes=sub_network,
                                            title=f"Enriched pathways of targets of {mirna}", msigdb=self.msigdb)
        image_file = html_file.split('/')[-1]
        html_image_path = f'{image_file}_PathwayEnriched.png'
        image_file = f'{self.directory}/{image_file}_PathwayEnriched.png'
        fig.savefig(image_file, dpi=300, bbox_inches='tight')

        with open(html_file, "r") as file:
            html_content = file.read()

        # Append the image tag to the HTML
        html_content += f"""
                    <h2>Enriched Pathways</h2>
                    <img src="{html_image_path}" alt="Enriched Pathways">
                    """

        # Save the modified HTML file
        with open(html_file, "w") as file:
            file.write(html_content)
        return genes
    def get_mirna_impact_in_cluster(self, mirnas):
        sub_network_cluster = []
        impacts = []
        for mirna in mirnas:
            sub_network, impact = self.get_mirna_influence(mirna)
            sub_network_cluster.extend(sub_network)
            impacts.append(impact)
        sub_network_cluster = list(set(sub_network_cluster))
        impacts = pd.DataFrame(impacts)
        genes = [gene for gene in sub_network_cluster if gene in impacts.columns]
        impacts = impacts[genes]
        impacts = impacts.reset_index()
        return impacts
    def create_report_clusters(self, save_directory=None):
        if save_directory is None:
            save_directory = self.directory
        for cluster, mirnas in self.clustered_mirnas.items():
            sub_network_cluster = []
            impacts = []
            for mirna in mirnas:
                sub_network, impact = self.get_mirna_influence(mirna)
                sub_network_cluster.extend(sub_network)
                impacts.append(impact)
            sub_network_cluster = list(set(sub_network_cluster))
            impacts = pd.DataFrame(impacts)
            genes = [gene for gene in sub_network_cluster if gene in impacts.columns]
            impacts = impacts[genes]
            impacts = impacts.reset_index()
            mi_plot.generate_html_network_report(network=self.mirna_network.network, selected_nodes=sub_network_cluster,
                                                 name=f"cluster_{cluster}", influences=impacts.round(3),
                                                 save_path=save_directory)
            # mi_plot.draw_network(G=my_network.network, node_list = sub_network_cluster, name= mirna+f"_cluster_{cluster}", save_path='mirna_scoring/sub_plots')
            html_file = f"{save_directory}/cluster_{cluster}.html"
            pathways_df = self.get_enriched_pathways(selected_genes=sub_network_cluster, fdr_threshold = 0.1)
            pathways_df.to_csv(f"{save_directory}/cluster_{cluster}_pathways.csv")
            pathways_df = pathways_df[pathways_df['FDR p-value']<0.05]
            pathways_df["collection"] = pathways_df["Term"].map(get_collection)
            for collection in pathways_df["collection"].unique():
                pathway_df = pathways_df[pathways_df['collection']==collection]
                if len(pathway_df)>0:
                    fig, ax = plot_ora_results(pathway_df, top_n=10, figsize=(12, 6), scale_odds_ratio=.5,
                             fontsize_title=12, fontsize_subtitle=12, fontsize_text=10,title=f"{collection} enriched pathways of targets of cluster {cluster}")
                    image_file = html_file.split('/')[-1]
                    html_image_path = f'{image_file}_PathwayEnriched_{collection}.png'
                    image_file = f'{save_directory}/{image_file}_PathwayEnriched_{collection}.png'
                    fig.savefig(image_file, dpi=300, bbox_inches='tight')

                    with open(html_file, "r") as file:
                        html_content = file.read()

                    # Append the image tag to the HTML
                    html_content += f"""
                        <h2>Enriched Pathways {collection}</h2>
                        <img src="{html_image_path}" alt="Enriched Pathways">
        
                        <h4> Filter to see in cytoscape: </h4>
                        """
                    with open(html_file, "w") as file:
                        file.write(html_content)
            html_content += jf.get_cytoscape_filter_from_list(sub_network_cluster)

            # Save the modified HTML file
            with open(html_file, "w") as file:
                file.write(html_content)

    def get_enriched_pathways(self, selected_genes, fdr_threshold = 0.1):
        enriched = dc.get_ora_df(
            df=selected_genes,
            net=self.msigdb,
            source='geneset',
            target='genesymbol'
        )
        pathway_df = enriched[enriched['FDR p-value'] < fdr_threshold]
        pathway_df.index = pathway_df["Term"]
        # pathway_df.set_index("Term", inplace=True)  # Set "Term" as index
        return pathway_df






