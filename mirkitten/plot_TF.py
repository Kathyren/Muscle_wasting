import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
from anndata import AnnData
import yaml

class plot_RF:
    def __init__(self, comparison, tf_df):
        """
        Initialize the plot_RF class with experiment and metadata file paths.
        """
        self.comparison = comparison
        self.tf_df = tf_df
    def plot_personalize_bars(self, title="", top=20, 
                            lFCs_thr=0.5, sig_thr=0.05,
                            sig = 'padj', title_fontsize=14, 
                            figsize=(5, 5), vertical=True, dpi=700, filename=None):
        # Generate the volcano plot
        fig = dc.plot_barplot(
            acts=self.tf_df,
            contrast=self.comparison,
            top=top,
            vertical=vertical,
            figsize=figsize,
            dpi=dpi,
            return_fig=True,
            save=filename
        )
        
        # Access the Axes object from the figure
        ax = fig.gca()
        
        # Add a title and remove the outer black line
        ax.set_title(title, fontsize=title_fontsize)
        for spine in ax.spines.values():
            spine.set_visible(False)  # Hide all outer spines
        
        return fig, ax
