import bisect
import numpy as np
from scipy import optimize
import scipy.stats as stats

from MethAnalysisHelpers.GenerateForm import GenerateForm

basic_func = GenerateForm().functional_form


def get_age_bins(ages, bin_number=20):
    binning = []
    interval = (max(ages) - min(ages)) / bin_number
    bin_ranges = [interval * x for x in range(bin_number)]
    for age in ages:
        binning.append(bisect.bisect_left(bin_ranges, age))
    return binning


def r2(x, y):
    # return r squared
    return stats.pearsonr(x, y)[0] ** 2


def mae(x, y):
    # return r squared
    return np.mean(abs(x - y))


def fit_trend(known: np.ndarray, predicted: np.ndarray, func=basic_func):
    popt, pcov = optimize.curve_fit(func, known, predicted)
    expected = func(known, *popt)
    return popt, pcov, expected


def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=1)
    phenotype_mean = np.mean(phenotype)

    # subtract means from observed values
    transformed_matrix = meth_matrix - matrix_means.reshape([-1, 1])
    transformed_phenotype = phenotype - phenotype_mean

    # calculate covariance
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))

    return covariance / (variance_meth * variance_phenotype)