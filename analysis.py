"""

OVERVIEW:  

Miscellanous analysis tools and scripts for general purpose data analysis.

"""

import warnings

def jsd(x,y): 
    # Jensen-shannon divergence
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log2(2*x/(x+y))
    d2 = y*np.log2(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)
    return d

def FDR_thresholds(q_threshold, N):
    # Benjamini - Hochberg criterion FDR thresholds
    # Use to reject H_0 in 1, ..., k for all p(i) <= (i/N)*(q/c(N))
    # where c(N) = sum_i 1/i
    #
    # Inputs: 
    #        q: q-value threshold (e.g. 0.05)
    #        N: number of tests
    # Returns:
    #        array of thresholds to be compared against ranked p-values
    cN = 0
    for i in range(1, N+1):
        cN += 1/float(i)
    FDR_tresholds = [(i/float(N))*(q_value/float(cN)) for i in range(1,N+1)]
    return FDR_thresholds


def corr_ratio(values, category_labels):
    # Computes correlation ratio for a given array of values and category_labels
    # Eqn: eta^2 = (sum_x [N_x * (mean(y_x) - mean(y))^2]) / (sum_x [sum_i [(y_xi - mean(y))^2]])
    # values: 1D array
    # category_labels: 1D array
    category_labels = np.array(category_labels)
    values = np.array(values)
    categories = np.unique(category_labels)
    data_dict = {}
    for catname in categories:
        indices = np.where(category_labels == catname)
        data_dict[catname] = values[indices]
    cat_means = {catname: np.mean(data_dict[catname]) for catname in categories}
    overall_mean = np.mean(values)
    nominator = np.sum([len(data_dict[catname])*(cat_means[catname]-overall_mean)**2 for catname in categories])
    denominator = np.sum([(val-overall_mean)**2 for catname in categories for val in data_dict[catname]])
    corr_ratio = nominator / denominator
    return corr_ratio
