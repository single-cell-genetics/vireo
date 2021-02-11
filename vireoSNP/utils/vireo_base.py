import numpy as np
from scipy.stats import entropy
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp, digamma, betaln, binom, gammaln


def get_binom_coeff(AD, DP, max_val=700, is_log=True):
    """Get the binomial coefficients
    """
    # Since binom can't give log value, the maximum value in 64bit is 
    # around e**700, close to binom(1000, 500)
    
    # print("Warning: this function is deprecated, please use logbincoeff.")
    idx = DP > 0
    _AD = AD[idx].astype(np.int64)
    _DP = DP[idx].astype(np.int64)
    
    binom_coeff = np.log(binom(_DP, _AD))
    binom_coeff[binom_coeff > max_val] = max_val
    binom_coeff = binom_coeff.astype(np.float32)
        
    return binom_coeff


def logbincoeff(n, k, is_sparse=False):
    """
    Ramanujan's approximation of log [n! / (k! (n-k)!)]
    This is mainly for convinience with pen. Please use betaln or gammaln
    """
    if is_sparse:
        RV_sparse = n.copy() * 0
        idx = (k > 0).multiply(k < n)
        n = np.array(n[idx]).reshape(-1)
        k = np.array(k[idx]).reshape(-1)
        
    RV = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    
    if is_sparse:
        RV_sparse[idx] += RV
        RV = RV_sparse
    return RV


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


def beta_entropy(X, X_prior=None, axis=None):
    """
    Get the entropy for beta distributions. If X_prior is not None, return the
    Kullback-Leibler divergence
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
    https://en.wikipedia.org/wiki/Beta_distribution#Quantities_of_information_(entropy)

    Parameters
    ----------
    X, X_prior: 
        numpy.array with shape: (N, 2)

    Example
    -------
    theta_shapes1 = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    theta_shapes2 = np.array([[364, 24197], [5886, 7475], [6075, 397]])
    beta_entropy(theta_shapes2)
    beta_entropy(theta_shapes2, theta_shapes1)
    """
    def _beta_cross_entropy(Xp, Xq):
        """return cross entropy -E_p[log q] for beta distribution
        For entropy, use as _beta_cross_entropy(X, X)
        """
        return (
            betaln(Xq[:, 0], Xq[:, 1]) - 
            (Xq[:, 0] - 1) * digamma(Xp[:, 0]) - 
            (Xq[:, 1] - 1) * digamma(Xp[:, 1]) +
            (Xq.sum(axis=1) - 2) * digamma(Xp.sum(axis=1))
            )

    # check shape
    if len(X.shape) == 1:
        if X.shape[0] == 2:
            X = X.reshape(-1, 2)
        else:
            print("Error: unsupported shape. Make sure it's (N, 2)")

    if X_prior is not None and len(X.shape) == 1:
        if X_prior.shape[0] == 2:
            X_prior = X_prior.reshape(-1, 2)
        else:
            print("Error: unsupported shape. Make sure it's (N, 2)")

    if X_prior is None:
        # entropy
        RV_mat = _beta_cross_entropy(X, X)
    else:
        # KL divergence
        RV_mat = _beta_cross_entropy(X, X_prior) - _beta_cross_entropy(X, X)

    return np.sum(RV_mat, axis=axis)


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


def optimal_match(X, Z, axis=1, return_delta=False):
    """
    Match Z to X by minimize the difference, 
    hence np.take(Z, idx1, aixs) is best aligned to np.take(X, idx0, aixs)
    
    Hungarian algorithm is used: 
    https://docs.scipy.org/doc/scipy-1.4.0/reference/generated/scipy.optimize.linear_sum_assignment.html
    """
    X_copy = X.copy()
    Z_copy = Z.copy()
    diff_mat = np.zeros((X.shape[axis], Z.shape[axis]))
    for i in range(X.shape[axis]):
        for j in range(Z.shape[axis]):
            diff_mat[i, j] = np.mean(np.abs(np.take(X_copy, i, axis=axis) - 
                                            np.take(Z_copy, j, axis=axis)))        
    idx0, idx1 = linear_sum_assignment(diff_mat)    
    if return_delta:
        return idx0, idx1, diff_mat
    else:
        return idx0, idx1


def greed_match(X, Z, axis=1):
    """
    This method has been dispatched, please use optimal_match!
    """
    print("This method has been dispatched, please use optimal_match!")
    return optimal_match(X, Z, axis=axis)[1]


def donor_select(GT_prob, ID_prob, n_donor, mode="distance"):
    """
    Select the donors from a set with extra donors.

    The GT_prior can have different number of donors from n_donor.
    
    mode="size": only keep the n_donor with largest number of cells
    mode="distance": only keep the n_donor with most different GT from each other
    """
    _donor_cnt = np.sum(ID_prob, axis=0)
    if mode == "size":
        _donor_idx = np.argsort(_donor_cnt)[::-1]
    else:
        _GT_diff = np.zeros((GT_prob.shape[1], GT_prob.shape[1]))
        for i in range(GT_prob.shape[1]):
            for j in range(GT_prob.shape[1]):
                _GT_diff[i, j] = np.mean(np.abs(GT_prob[:, i, :] - 
                                                GT_prob[:, j, :]))

        _donor_idx = [np.argmax(_donor_cnt)]
        _donor_left = np.delete(np.arange(GT_prob.shape[1]), _donor_idx)
        _GT_diff = np.delete(_GT_diff, _donor_idx, axis=1)
        while len(_donor_idx) < _GT_diff.shape[0]:
            # _idx = np.argmax(np.sum(_GT_diff[_donor_idx, :], axis=0))
            _idx = np.argmax(np.min(_GT_diff[_donor_idx, :], axis=0))
            _donor_idx.append(_donor_left[_idx])
            _donor_left = np.delete(_donor_left, _idx)
            _GT_diff = np.delete(_GT_diff, _idx, axis=1)

    print("[vireo] donor size with searching extra %d donors:" 
              %(GT_prob.shape[1] - n_donor))
    print("\t".join(["donor%d" %x for x in _donor_idx]))
    print("\t".join(["%.0f" %_donor_cnt[x] for x in _donor_idx]))

    ID_prob_out = ID_prob[:, _donor_idx[:n_donor]]
    ID_prob_out[ID_prob_out < 10**-10] = 10**-10

    return ID_prob_out
