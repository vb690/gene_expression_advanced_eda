from copy import deepcopy

import numpy as np


class CatEncoderDecoder:
    """
    """
    def __init__(self, categories):
        """
        """
        self.cat_to_num = {cat: num for num, cat in enumerate(categories)}
        self.num_to_cat = {num: cat for cat, num in self.cat_to_num.items()}

        self.cat_to_oh = {}
        for num, cat in enumerate(categories):

            oh = np.zeros(shape=(len(categories)))
            oh[num] = 1
            self.cat_to_oh[cat] = oh

    def encode(self, to_encode, num=True):
        """
        """
        if num:
            encoded = np.array([self.cat_to_num[cat] for cat in to_encode])
        else:
            encoded = np.array([self.cat_to_oh[cat] for cat in to_encode])
        return encoded

    def decode(self, to_decode, num=True):
        """
        """
        if num:
            decoded = np.array([self.num_to_cat[num] for num in to_decode])
        else:
            decoded = np.array([np.argmax(oh) for oh in to_decode])
        return decoded


def cohen_d(t, n):
    """Utility function for computing Cohen's D given t statistics and n
    """
    d = (2*t) / ((n-1) ** 0.5)
    return d


def log2_fold_change(term_1, term_2):
    """Utility fucntion for computing log2 fold change
    """
    log_2_term_1 = np.log2(
        term_1,
        where=term_1 > 0
    )
    log_2_term_2 = np.log2(
        term_2,
        where=term_2 > 0
    )
    term_1_mean = np.nanmean(log_2_term_1, axis=0)
    term_2_mean = np.nanmean(log_2_term_2, axis=0)
    fold_change = term_1_mean - term_2_mean
    return fold_change


def top_k_variance(X, names=None, k=3000, no_variance_filter=True):
    """Utility function for retaining top k features with the highest variance.
    """
    var = np.var(X, axis=0)
    if no_variance_filter:
        mask = var != 0
        X = X[:, mask]
        var = np.var(X, axis=0)
        names = names[mask]
    sorted_var = np.sort(var)[::-1]
    threshold = sorted_var[k - 1]

    mask = var >= threshold
    if names is not None:
        return X[:, mask], names[mask]
    else:
        return X[:, mask]


def introduce_outliers(expression_values,  gene_name, outliers_perc=0.1,
                       verbose=True):
    """Utility function for introducing artificial outliers
    in the samples.

    An outlier will be sampled from a normal distribution having mean
    as the mean of a specific gene expression values + or - 2.5 times its
    standard deviation.

    Args:
        - expression_values (iterable): expression values for a specific gene
        - gene_name (str): name of the gene
        - outliers_perc (float): percentage of outlying samples introduced
        - verbose (bool): bolean, control the verbosity of the function

    Returns:
        - new_name (str): new name indicating outliers have been introduced
        - new_expression_values (iterable): expression_values with outliers
    """
    expression_values = deepcopy(expression_values)
    exp_mean = np.nanmean(expression_values)
    exp_std = np.nanstd(expression_values)
    new_name = 'outlying_' + gene_name

    samples = [sample for sample in range(len(expression_values))]
    samples = np.random.choice(
        a=samples,
        size=int(len(samples) * outliers_perc),
        replace=False
    )

    for sample_outlier in samples:

        if verbose:
            print('')
            print(f'Adding outliers to sample {sample_outlier}:')

        if np.random.choice(['plus', 'minus']) == 'minus':
            # we can't have negative gene expression
            outlier_value = np.random.normal(
                loc=exp_mean - (2.5 * exp_std),
                scale=exp_std
            )
        else:
            outlier_value = np.random.normal(
                loc=exp_mean + (2.5 * exp_std),
                scale=exp_std
            )

        expression_values[sample_outlier] = outlier_value

    return new_name, expression_values
