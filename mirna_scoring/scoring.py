
import mirna_scoring.mirna_impact as mis
import mirna_scoring.mirna_influence_plots as mi_plot
from network.network_processing import load_graph


class scoring:
    def __init__(self, network_path, steps, keywords, samples, save_name, safe_path='results/'):
        self.network_path = network_path
        self.steps = steps
        self.keywords = keywords
        self.samples = samples
        self.save_name = save_name
        self.network = load_graph(f"network/Networks_pkl/complete_n_tf_mirnas__UseCase5_cutoff_0.95.pkl")
    def score_mirnas(self):
        """
        Perform microRNA scoring on the network.
        """
        mirna_network=mis.mirna_network(self.network)
        evaluation = mis.mirna_evaluation(mirna_network=mirna_network, session_name=self.save_name)
        evaluation_df = evaluation.score(steps=self.steps, sample_size=self.sample_size, pathway_keywords=self.pathway_keywords)
        return evaluation_df
    def save_results(self, evaluation_df):
        """
        Save the evaluation results to a CSV file.
        """
        evaluation_df.to_csv(f"{self.save_name}.csv", index=False)
        print(f"Results saved to {self.save_name}.csv")
