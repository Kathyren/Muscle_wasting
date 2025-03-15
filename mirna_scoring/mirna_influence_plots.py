import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import cityblock

def convert_to_int_list(lst):
    """
    This function is only for specific formating,
    the output of impact is a list as [1,1,-1,1] etc. This value can be read as a string
    when read the csv files, and this is just to convert the string to a list
    :param lst:
    :return:
    """
    if isinstance(lst, list):
        return [int(x) for x in lst]
    elif isinstance(lst, str):
        string = lst.replace('[', '').replace(']', '')
        string = string.split(',')
        if len(string) < 1:
            return []
        return [int(x) if x != '' else 0 for x in string]
    else:
        return [lst]



def plot_dotplot(df, gene_scale=0.3, mirna_scale=1.5):
    """
    This take the dataframe with the lists of [1,1,-1] type, and plots a dotplot
    with the size of the dot the len of the list and the color the sum
    :param df:
    :return:
    """
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
    n_mid = abs(min(plot_data['Influence'])) / (abs(min(plot_data['Influence'])) + max(plot_data['Influence']))
    nodes = [0, n_mid, 1]
    cmap = LinearSegmentedColormap.from_list('custom', list(zip(nodes, colors)))

    # Plot using seaborn
    n_mirna = len(size_data.columns) + 1
    n_genes = len(size_data) + 1
    print('genes', n_genes - 1)
    print('mirnas', n_mirna - 1)
    height_plot = min(655, int(n_genes * gene_scale))

    height_plot = max(10, height_plot)
    sns.set(rc={'axes.facecolor': 'lightgray'})
    x = n_mirna * mirna_scale
    y = height_plot
    fig, ax = plt.subplots(figsize=(x, y))
    ax.grid(False)

    # plt.figure(figsize=(10, height_plot))
    scatter_plot = sns.scatterplot(data=plot_data, x='miRNA', y='Gene', size='Paths', hue='Influence', palette=cmap,
                                   sizes=(x, y))
    scatter_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    if n_genes > 2:
        # Get the first two and last y-tick positions.
        miny, nexty, *_, maxy = ax.get_yticks()

        # Compute half the y-tick interval (for example).
        eps = (nexty - miny) / 2  # <-- Your choice.
        plt.xticks(rotation=90)
        # Adjust the limits.
        ax.set_ylim(maxy + eps, miny - eps)
    # plt.tight_layout()
    plt.show()
    return fig,ax


def plot_pathway_frequency(frequency_pathways):

    x = frequency_pathways.keys()
    y = frequency_pathways.values()
    plt.figure(figsize=(20,6))
    plt.xticks(rotation=90)

    plt.plot(x, y)
    plt.show()

def plot_pathways_keyword_heatmap(mir_pathway_influence_df):
    mir_pathway_influence_df['participation'] = mir_pathway_influence_df.drop(
        columns=["Different_pathways", "Total"]).sum(axis=1)
    mir_pathway_influence_df_n0 = mir_pathway_influence_df[mir_pathway_influence_df['participation'] > 0]
    mir_pathway_influence_df_n0 = mir_pathway_influence_df_n0.drop(
        columns=["Different_pathways", "Total", "participation"])
    mir_pathway_influence_df_n0 = mir_pathway_influence_df_n0.loc[:, (mir_pathway_influence_df_n0 != 0).any(axis=0)]
    sns.heatmap(mir_pathway_influence_df_n0, cmap="YlOrBr", annot=True)


def plot_mirnas_similarirty(dist_df):
    sorted_index = dist_df.mean().sort_values().index  # Sorting by mean distance

    # Apply sorting to the DataFrame
    sorted_df = dist_df.loc[sorted_index, sorted_index]

    # Plot the sorted distance matrix as a heatmap
    n_mirnas = len(sorted_df)
    size_l = n_mirnas * 0.3
    plt.figure(figsize=(size_l, size_l))
    sns.heatmap(sorted_df, annot=False, cmap='coolwarm')  # , linewidths=0.5)#, linecolor='black')
    plt.title('Sorted Distance Matrix Heatmap')
    plt.show()