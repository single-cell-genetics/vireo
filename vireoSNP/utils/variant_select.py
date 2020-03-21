import numpy as np
from scipy.stats import entropy

def barcode_entropy(X, y=None):
    """
    entropy for categorical barcodes
    """
    if y is None:
        Z_str = [str(x) for x in X]
    elif len(X) == len(y):
        Z_str = [str(X[i]) + str(y[i]) for i in range(len(X))]
    else:
        print("Error: X and y have different length in barcode_entropy.")
        return None, None
    
    Z_val, Z_cnt = np.unique(Z_str, return_counts=True)
    
    return entropy(Z_cnt / np.sum(Z_cnt), base=2), Z_str


def variant_select(GT, var_count=None, rand_seed=0):
    """
    Selection of a set of discriminatory variants by prioritise variants on 
    information gain.

    GT: (n_var * n_donor)
        a matrix with categorical values
    var_count: (n_var, )
        the counts for each variant
    """
    np.random.seed(rand_seed)
    
    K = GT.shape[1]
    entropy_now = 0
    variant_set = []
    barcode_set = ["#"] * K

    entropy_all = np.zeros(GT.shape[0])
    barcode_all = [barcode_set] * GT.shape[0]
    while True:
        for i in range(GT.shape[0]):
            _entropy, _barcode = barcode_entropy(barcode_set, GT[i, :])
            entropy_all[i], barcode_all[i] = _entropy, _barcode
        if np.max(entropy_all) == entropy_now:
            break
        
        idx = np.where(np.max(entropy_all) == entropy_all)[0]
        if var_count is not None:
            # idx = idx[np.argsort(var_count[idx])[::-1]]
            idx = idx[var_count[idx] >= np.median(var_count[idx])]
        print("Randomly select 1 more variants out %d" %len(idx))
        idx_use = idx[np.random.randint(len(idx))]
        
        variant_set.append(idx_use)
        barcode_set = barcode_all[idx_use]
        entropy_now = entropy_all[idx_use]
        
    if entropy_now < np.log2(K):
        print("Warning: variant_select can't distinguish all samples.")

    return entropy_now, barcode_set, variant_set
