import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import seaborn as sns

def plot_dotplot(df):
    def convert_to_int_list(lst):
        return [int(x) for x in lst]

    # Create new DataFrame for plotting
    plot_data = pd.DataFrame(index=df.index, columns=df.columns)
    size_data = pd.DataFrame(index=df.index, columns=df.columns)
    color_data = pd.DataFrame(index=df.index, columns=df.columns)

    for col in df.columns:
        for idx in df.index:
            int_list = convert_to_int_list(df.at[idx, col])
            size_data.at[idx, col] = len(int_list)
            color_data.at[idx, col] = sum(int_list)

    # Convert to long format for seaborn
    plot_data = pd.DataFrame({
        'Gene': np.repeat(df.index, df.shape[1]),
        'miRNA': np.tile(df.columns, df.shape[0]),
        'Paths': size_data.values.flatten(),
        'Influence': color_data.values.flatten()
    })

    # Define custom colormap
    cmap = LinearSegmentedColormap.from_list('custom', ['red', 'white', 'blue'])

    # Filter out rows where the size is 0
    plot_data = plot_data[plot_data['Paths'] > 0]

    colors = ['red', 'white', 'blue']
    n_mid = abs(min(plot_data['Influence'])) / (abs(min(plot_data['Influence'])) + abs(max(plot_data['Influence'])))
    nodes = [0, n_mid, 1]
    cmap = LinearSegmentedColormap.from_list('custom', list(zip(nodes, colors)))

    # Plot using seaborn
    n_genes = len(plot_data)
    print(n_genes)
    height_plot = min(655, int(n_genes * .21))
    sns.set(rc={'axes.facecolor': 'lightgray'})
    fig, ax = plt.subplots(figsize=(15, height_plot))
    ax.grid(False)
    # plt.figure(figsize=(10, height_plot))
    scatter_plot = sns.scatterplot(data=plot_data, x='miRNA', y='Gene', size='Paths', hue='Influence', palette=cmap,
                                   sizes=(20, 200))
    scatter_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

    # Get the first two and last y-tick positions.
    miny, nexty, *_, maxy = ax.get_yticks()

    # Compute half the y-tick interval (for example).
    eps = (nexty - miny) / 2  # <-- Your choice.

    # Adjust the limits.
    ax.set_ylim(maxy + eps, miny - eps)
    # plt.tight_layout()
    plt.show()


def plot_pathway_frequency(frequency_pathways):

    x = frequency_pathways.keys()
    y = frequency_pathways.values()
    plt.figure(figsize=(20,6))
    plt.xticks(rotation=90)

    plt.plot(x, y)
    plt.show()