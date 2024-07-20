import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize

from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering


# Function to ensure all elements are integers
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
        return [int(x) for x in string]
    else:
        return [lst]


def get_impact_data(df):
    """
    From the list of impact, it just sums the numbers and return a dataframe where each value is an integer
    :param df:
    :return:
    """
    color_data = pd.DataFrame(index=df.index, columns=df.columns)

    for col in df.columns:
        for idx in df.index:
            int_list = convert_to_int_list(df.at[idx, col])
            color_data.at[idx, col] = sum(int_list)

    return color_data


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
    fig, ax = plt.subplots(figsize=(x,y ))
    ax.grid(False)

    # plt.figure(figsize=(10, height_plot))
    scatter_plot = sns.scatterplot(data=plot_data, x='miRNA', y='Gene', size='Paths', hue='Influence', palette=cmap,
                                   sizes=(x, y))
    scatter_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

    # Get the first two and last y-tick positions.
    miny, nexty, *_, maxy = ax.get_yticks()

    # Compute half the y-tick interval (for example).
    eps = (nexty - miny) / 2  # <-- Your choice.
    plt.xticks(rotation=90)
    # Adjust the limits.
    ax.set_ylim(maxy + eps, miny - eps)
    # plt.tight_layout()
    plt.show()


def order_length_mirna(df, interest_mirna: str):
    """
    With the original df with the lists, get the length of the list, and order the dataframe given one specific mirna
    The df has mirnas on the columns and genes in the rows
    :param df:
    :param interest_mirna:
    :return:
    """
    df_sorted = df.copy()
    df_sorted['length'] = df_sorted[interest_mirna].apply(lambda x: len(x) if isinstance(x, list) else 0)
    df_sorted = df_sorted.sort_values(by='length', ascending=False)
    df_sorted = df_sorted.drop(columns=['length'])
    return df_sorted


def order_impact_mirna(df, interest_mirna: str):
    """
     With the original df with the lists, get the sum of the list, and order the dataframe given one specific mirna
    The df has mirnas on the columns and genes in the rows
    :param df:
    :param interest_mirna:
    :return:
    """
    df_sorted = df.copy()
    df_sorted['impact'] = df_sorted[interest_mirna].apply(lambda x: sum(x) if isinstance(x, list) else 0)
    df_sorted = df_sorted.sort_values(by='impact', ascending=False)
    df_sorted = df_sorted.drop(columns=['impact'])
    return df_sorted


# Function to calculate the sum of lengths of lists in a row
def sum_lengths(row):
    """
    Sum the values of the list
    :param row:
    :return:
    """
    return sum(len(lst) for lst in row)


def get_minra_influence(test_mir, df):
    """

    :param test_mir:
    :param df:
    :return:
    """
    test_mir_df = df.T[test_mir]
    test_mir_df = pd.DataFrame(test_mir_df)
    test_mir_df['length'] = test_mir_df[test_mir].apply(lambda x: len(x) if isinstance(x, list) else 0)
    test_mir_df = test_mir_df.sort_values(by='length', ascending=False)
    test_mir_df = test_mir_df.drop(columns=['length'])
    return test_mir_df


def get_mirnas_similar_impact(df):
    """
    This calculates the manhattan distance of the microRNAs impact on the genes. Getting the
    similarity of the mirnas given thei impact
    :param df:
    :return:
    """
    df_numeric = df.apply(pd.to_numeric, errors='coerce')

    # Calculate the pairwise distances between columns
    dist_matrix = pdist(df_numeric.T, metric='cityblock')
    dist_matrix_square = squareform(dist_matrix)

    # Create a DataFrame for the distance matrix
    dist_df = pd.DataFrame(dist_matrix_square, index=df_numeric.columns, columns=df_numeric.columns)
    sorted_index = dist_df.mean().sort_values().index  # Sorting by mean distance

    return dist_df.loc[sorted_index, sorted_index]


def plot_similarity_heatmap(sorted_df):
    """
    Plot the sorted distance matrix as a heatmap

    :param sorted_df: the distance matrix of the microRNA and their distance
    :return:
    """
    #
    n_mirnas = len(sorted_df)
    size_l = n_mirnas * 0.3
    plt.figure(figsize=(size_l, size_l))
    sns.heatmap(sorted_df, annot=False, cmap='coolwarm')  # , linewidths=0.5)#, linecolor='black')
    plt.title('Sorted Distance Matrix Heatmap')
    plt.show()


def cluster_mirnas(dist_matrix_square, n_clusters=10):
    """
    This function takes the distance matrix of the micrornas and cluster them by similarity
    return the list of the micrornas and their cluster.
    :param dist_matrix_square:
    :param n_clusters:
    :return:
    """
    # Convert distance matrix to similarity matrix (affinity matrix)
    affinity_matrix = np.exp(-dist_matrix_square / dist_matrix_square.std())  # Using Gaussian kernel for similarity

    # Apply Spectral Clustering
    # Number of clusters
    sc = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=42)
    cluster_labels = sc.fit_predict(affinity_matrix)

    # Print cluster labels
    print("Cluster labels:\n", cluster_labels)
    mirnas = dist_matrix_square.index

    mirna_clusters = pd.DataFrame(cluster_labels, index=mirnas, columns=['Cluster'])
    return mirna_clusters


def get_curve_type(de_df):
    """
    Given the DE log fold change od the genes, it calculates the type of changes the gene has on life.
    A means the gene didn't change
    B that the gene decreased with age
    C that the gene increased with age

    It evaluates for young to middle age, and middle age to old. And finally, it adds -A, -B, -C depending on the
    change fotm young to old

    For example, if the gene is upregulated in young vs middle age, down regulated in middle age vs old, this means
    BC

    and from there, if there is differentially express change yong vs old,
    it can be BC-B if there is more in young vs old
    :param de_df:
    :return:
    """
    curves = []
    if 'ym' in de_df.columns and \
            'mo' in de_df.columns and \
            'yo' in de_df.columns:

        for gene, row in de_df.iterrows():
            curve_type = ['X', 'X']
            if pd.isna(row['ym']):
                curve_type[0] = 'A'
            elif row['ym'] > 0:
                curve_type[0] = 'B'
            elif row['ym'] < 0:
                curve_type[0] = 'C'
            if pd.isna(row['mo']):
                curve_type[1] = 'A'
            elif row['mo'] > 0:
                curve_type[1] = 'B'
            elif row['mo'] < 0:
                curve_type[1] = 'C'
            #if curve_type == ['C', 'B'] or curve_type == ['B', 'C']:
            if pd.isna(row['yo']):
                curve_type.append('')
            elif row['yo'] > 0:
                curve_type.append('-B')
            else:
                curve_type.append('-C')
            curve_type = ''.join(curve_type)
            curves.append(curve_type)
    return curves


def calculate_measurements(int_influence_df):
    influence_weight = int_influence_df.copy()
    influence_weigh_ym = int_influence_df.copy()
    influence_weigh_mo = int_influence_df.copy()
    influence_weigh_yo = int_influence_df.copy()
    influence_quantity = int_influence_df.copy()
    for col in int_influence_df.columns:
        if col not in ['yo', 'ym', 'mo']:
            yo = int_influence_df[col] * int_influence_df['yo']
            ym = int_influence_df[col] * int_influence_df['ym']
            mo = int_influence_df[col] * int_influence_df['mo']
            influence_weigh_ym[col] = ym
            influence_weigh_mo[col] = mo
            influence_weigh_yo[col] = yo
            yo = (int_influence_df[col] * int_influence_df['yo'] > 0).astype(int)
            ym = (int_influence_df[col] * int_influence_df['ym'] > 0).astype(int)
            mo = (int_influence_df[col] * int_influence_df['mo'] > 0).astype(int)
            influence_quantity[col] = yo + ym + mo

    influence_weight.drop(columns=['yo', 'ym', 'mo'], inplace=True)
    influence_weigh_ym.drop(columns=['yo', 'ym', 'mo'], inplace=True)
    influence_weigh_mo.drop(columns=['yo', 'ym', 'mo'], inplace=True)
    influence_weigh_yo.drop(columns=['yo', 'ym', 'mo'], inplace=True)

    influence_quantity.drop(columns=['yo', 'ym', 'mo'], inplace=True)
    measurements = [influence_weight, influence_weigh_ym, influence_weigh_mo, influence_weigh_yo, influence_quantity]
    return measurements
