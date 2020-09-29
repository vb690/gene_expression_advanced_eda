from copy import deepcopy

from itertools import combinations

import math

import numpy as np

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests as mt

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


def visualize_PCA_summary(pca_fitter, features_name, title, k=40,
                          figsize=(10, 20),  **kwargs):
    """Utility function for visualizing the top and bottom k features
    in terms of their association with the first 2 principal components.

    Args:
        - pca_fitter (PCA object): a fitted PCA object
        - features_name(iterable): names of the features passed to PCA
        - title(string): general title passed to the plot
        - k(int): top and bottom k features
                 in terms of their association with the first 2 principal
                 components

    Returns:
        - None
    """
    df = pd.DataFrame(
        columns=['component_1', 'component_2', 'features_name']
    )
    df['component_1'] = pca_fitter.components_[0, :]
    df['component_2'] = pca_fitter.components_[1, :]
    df['features_name'] = features_name

    fig, axs = plt.subplots(1, 2, figsize=figsize)
    fig.suptitle(title)
    for ax, n_component in zip(axs.flatten(), [1, 2]):

        component_name = f'component_{n_component}'
        # Visualize weights for component 1
        df = df.sort_values(component_name)
        sns.barplot(
            y='features_name',
            x=component_name,
            data=df.iloc[np.r_[0:k, -k:0]],
            ax=ax,
            **kwargs
        )
        ax.set_xlabel(f'{component_name} Features Weights')

    plt.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.show()

    return None


def visualize_dim_reduction(reduction, title, outliers_loc=None, labels=None,
                            figsize=(10, 10), save_dir=None, **kwargs):
    """Utility function for visualizing the data in a lower dimensional space.

    No matter the number of components chosen
    the fucntion will plot only the first 2.

    Args:
        - reduction(numpy array): result of dimensionality reduction.
        - title(string): title for the plot
        - outliers_loc(iterable): index of outlying samples
        - labels(iterable): labels associated to each sample
        - **kwargs: keyword arguments passed to plt.scatter()

    Returns:
        - None
    """
    plt.figure(figsize=figsize)
    # if we have labels
    if labels is not None:
        unique_labels = np.unique(labels).flatten()

        for unique_label in unique_labels:

            indices = np.argwhere(labels == unique_label).flatten()
            plt.scatter(
                reduction[indices, 0],
                reduction[indices, 1],
                label=unique_label,
                ** kwargs
            )
    else:
        plt.scatter(
            reduction[:, 0],
            reduction[:, 1],
            ** kwargs
        )
    # if we know where the outliers are
    if outliers_loc is not None:

        for loc in outliers_loc:

            plt.scatter(
                reduction[loc, 0],
                reduction[loc, 1],
                c='r',
                ** kwargs
            )
            plt.annotate(
                loc,
                (reduction[loc, 0], reduction[loc, 1])
            )

    plt.xlabel(f'Component 1')
    plt.ylabel(f'Component 2')
    plt.title(title)
    plt.legend()

    plt.tight_layout()

    if save_dir is not None:
        plt.savefig(
            f'{save_dir}\\{title}.png'
        )
        plt.close()
    else:
        plt.show()
        plt.close()

    return None


def visualize_expression(clusters_expressions, figsize=(10, 25),
                         median=True, **kwargs):
    """

    Args:
        - clusters_expressions(dictionary): keys are cluster,
                                           values are arrays with the genes
                                           expression for the specific cluster.
        - figsize:
        - median:
        - **kwargs:

    Returns
        - None
    """
    if median:
        clusters_expressions = {
            cluster: np.nanmedian(expressions, axis=0) for
            cluster, expressions in clusters_expressions.items()
        }
    fig, axs = plt.subplots(
        len(clusters_expressions),
        1,
        figsize=figsize,
        sharex=True
    )
    clusters = list(clusters_expressions.keys())
    for cluster, ax in zip(clusters, axs.flatten()):

        data = clusters_expressions[cluster]
        if median:
            data = data.reshape(1, -1)
        mask = None
        title = f'Genes Expression in Cluster {cluster}'

        sns.heatmap(
            data,
            cmap='coolwarm',
            ax=ax,
            mask=mask,
            **kwargs
        )
        ax.set_title(title)
        ax.set_xlabel('Genes')
        ax.set_xticks([])

    plt.show()


def visualize_expression_comp(clusters_expressions, figsize=(10, 40),
                              p_thresh=0.01, fc_thresh=1.5, **kwargs):
    """
    """
    clusters_expressions = deepcopy(clusters_expressions)
    comparisons = combinations(
        iterable=list(clusters_expressions.keys()),
        r=2
    )
    comparisons = [element for element in comparisons]
    fig, axs = plt.subplots(
        len(comparisons),
        1,
        figsize=figsize,
        sharex=True
    )
    for comparison, ax in zip(comparisons, axs.flatten()):

        term_1 = comparison[0]
        term_2 = comparison[1]

        log_2_term_1 = np.log2(
            clusters_expressions[term_1],
            where=clusters_expressions[term_1] > 0
        )
        log_2_term_2 = np.log2(
            clusters_expressions[term_2],
            where=clusters_expressions[term_2] > 0
        )

        t, p = ttest_ind(
            log_2_term_1,
            log_2_term_2,
            axis=0
        )
        p[np.isnan(p)] = 1
        reject, corrected_p, sidak, bonf = mt(p, alpha=p_thresh)

        difference = np.nanmean(log_2_term_1, axis=0) \
            - np.nanmean(log_2_term_2, axis=0)

        colors = []
        for p_val, diff in zip(corrected_p, difference):
            if p_val > p_thresh or -fc_thresh < diff < fc_thresh:
                colors.append('#DCDCDC')
            elif p_val < p_thresh and diff > fc_thresh:
                colors.append('#FF0000')
            elif p_val < p_thresh and diff < -fc_thresh:
                colors.append('#0000FF')

        ax.scatter(
            difference,
            [- math.log10(p_val + 1e-15) for p_val in corrected_p],
            c=colors,
            **kwargs
        )
        ax.axhline(-math.log10(p_thresh), c='g', linestyle='--', alpha=0.5)
        ax.axvline(-fc_thresh, c='r', linestyle='--', alpha=0.5)
        ax.axvline(fc_thresh, c='r', linestyle='--', alpha=0.5)

        ax.set_title(f'T-test cluster {term_1} against cluster {term_2}')
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-log10(corrected p-vals)')

        # ax.set_xlim((-4, 4))

    plt.show()
