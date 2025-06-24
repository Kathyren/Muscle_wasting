
import mirna_scoring.mirna_impact as mis
import mirna_scoring.mirna_influence_plots as mi_plot
from network.network_processing import load_graph
import cytoscape as ct
import network.network_processing as ntp

class scoring:
    def __init__(self, network_path, steps, keywords, samples, save_name, safe_path='results/'):
        self.network_path = network_path
        # network path is a pkl file, so we need to load it
        if not network_path.endswith('.pkl'):
            raise ValueError("Network path must be a .pkl file")
        if not network_path:
            raise ValueError("Network path cannot be empty")
        self.steps = steps
        self.keywords = keywords
        self.samples = samples
        self.save_name = save_name
        self.network = load_graph(network_path)


    
    def score_mirnas(self):
        """
        Perform microRNA scoring on the network.
        """
        mirna_network=mis.mirna_network(self.network)
        evaluation = mis.mirna_evaluation(mirna_network=mirna_network, session_name=self.save_name)
        evaluation_df = evaluation.score(steps=self.steps, sample_size=self.samples, pathway_keywords=self.keywords)
        return evaluation_df
    def save_results(self, evaluation_df):
        """
        Save the evaluation results to a CSV file.
        """
        evaluation_df.to_csv(f"{self.save_name}.csv", index=True)
        print(f"Results saved to {self.save_name}.csv")
