import numpy as np
from scipy.stats import entropy
from scipy.special import logsumexp, digamma, betaln, binom


def get_binom_coeff(AD, DP, max_val=700, is_log=True):
    """Get the binomial coefficients
    """
    # Since binom can't give log value, the maximum value in 64bit is 
    # around e**700, close to binom(1000, 500)
    idx = DP > 0
    _AD = AD[idx].astype(np.int64)
    _DP = DP[idx].astype(np.int64)
    
    binom_coeff = np.log(binom(_DP, _AD))
    binom_coeff[binom_coeff > max_val] = max_val
    binom_coeff = binom_coeff.astype(np.float32)
        
    return binom_coeff


def normalize(X, axis=-1):
    """
    Normalization of tensor with sum to 1.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    tensor_normalize(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    X_sum = np.sum(X, axis=axis, keepdims=True)
    return X / X_sum

def tensor_normalize(X, axis=1):
    return normalize(X, axis)


def loglik_amplify(X, axis=-1):
    """
    Amplify the log likelihood matrix by subtract the maximum.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    X_max = np.max(X, axis=axis, keepdims=True)
    return X - X_max


def beta_entropy(X, X_prior=None):
    """
    Get the entropy for beta distributions. If X_prior is not None, return the
    Kullback-Leibler divergence
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html

    Example
    -------
    theta_shapes1 = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    theta_shapes2 = np.array([[364, 24197], [5886, 7475], [6075, 397]])
    beta_entropy(theta_shapes2)
    beta_entropy(theta_shapes2, theta_shapes1)
    """
    RV1 = 0
    if X_prior is None:
        X_prior = X.copy()
    else:
        RV1 = (- betaln(X[:, 0], X[:, 1]) + 
               (X[:, 0] - 1) * digamma(X[:, 0]) + 
               (X[:, 1] - 1) * digamma(X[:, 1]) - 
               (X.sum(axis=1) - 2) * digamma(X.sum(axis=1)))

    RV2 = (- betaln(X_prior[:, 0], X_prior[:, 1]) + 
           (X_prior[:, 0] - 1) * digamma(X[:, 0]) + 
           (X_prior[:, 1] - 1) * digamma(X[:, 1]) - 
           (X_prior.sum(axis=1) - 2) * digamma(X.sum(axis=1)))

    return np.sum(RV1) - np.sum(RV2)


def match(ref_ids, new_ids, uniq_ref_only=True):
    """
    Mapping new_ids to ref_ids. ref_ids can have repeated values, but new_ids 
    can only have unique ids or values. Therefore, new_ids[RT_idx] will be 
    the same as ref_ids. Note, 
    
    Parameters
    ----------
    ref_ids : array_like or list
        ids for reference with type of int, float, or string
    new_ids : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of ref_ids
        The index for new_ids mapped to ref_ids. If an id in ref_ids does not 
        exist in new_ids, then return a None for that id. 
    Examples
    --------
    >>> x1 = [5, 9, 1]
    >>> x2 = [1, 2, 5, 7, 9]
    >>> match(x1, x2)
    array([2, 4, 0])
    >>> match(x2, x1)
    array([2, None, 0, None, 1], dtype=object)
    >>> RT_idx = match(x2, x1)
    >>> idx1 = numpy.where(RT_idx != None)[0]
    >>> idx1
    array([0, 2, 4])
    >>> idx2 = RT_idx[idx1].astype(int)
    >>> idx2
    array([2, 0, 1])
    """
    idx1 = np.argsort(ref_ids)
    idx2 = np.argsort(new_ids)
    RT_idx1, RT_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or ref_ids[idx1[i]] < new_ids[idx2[j]]:
            RT_idx1.append(idx1[i])
            RT_idx2.append(None)
            i += 1
        elif ref_ids[idx1[i]] == new_ids[idx2[j]]:
            RT_idx1.append(idx1[i])
            RT_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: j += 1
        elif ref_ids[idx1[i]] > new_ids[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RT_idx1)
    RT_idx = np.array(RT_idx2)[origin_idx]
    return RT_idx
