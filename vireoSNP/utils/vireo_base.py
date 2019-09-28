import numpy as np
from scipy.stats import entropy
from scipy.special import logsumexp, digamma, betaln


def tensor_normalize(X, axis=1):
    """
    Normalization of tensor with sum to 1.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    tensor_normalize(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    X_sum = np.sum(X, axis=axis).reshape(shape2)
    return X / X_sum


def loglik_amplify(X, axis=1):
    """
    Amplify the log likelihood matrix by subtract the maximum.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    X_max = np.max(X, axis=axis).reshape(shape2)
    return X - X_max


def get_theta_shapes(AD, DP, ID_prob, GT_prob, theta_prior):
    """
    """
    S1_gt = AD * ID_prob
    SS_gt = DP * ID_prob
    S2_gt = SS_gt - S1_gt
    
    theta_shapes = theta_prior.copy()
    for ig in range(theta_shapes.shape[0]):
        theta_shapes[ig, 0] += np.sum(S1_gt * GT_prob[:, ig, :])
        theta_shapes[ig, 1] += np.sum(S2_gt * GT_prob[:, ig, :])
    return theta_shapes

def get_ID_prob(AD, DP, GT_prob, theta_shapes, Psi=None):
    """
    """
    if Psi is None:
        Psi = np.ones(GT_prob.shape[2]) / GT_prob.shape[2]
        
    logLik_ID = np.zeros((AD.shape[1], GT_prob.shape[2]))
    for ig in range(GT_prob.shape[1]):
        S1 = AD.transpose() * GT_prob[:, ig, :]
        SS = DP.transpose() * GT_prob[:, ig, :]
        S2 = SS - S1
        logLik_ID += (S1 * digamma(theta_shapes[ig, 0]) + 
                      S2 * digamma(theta_shapes[ig, 1]) - 
                      SS * digamma(np.sum(theta_shapes[ig, :])))
    
    Psi_norm = np.log(Psi / np.sum(Psi))
    ID_prob = np.exp(loglik_amplify(logLik_ID + Psi_norm, axis=1))
    ID_prob = tensor_normalize(ID_prob, axis=1)
    
    return ID_prob, logLik_ID
    

def get_GT_prob(AD, DP, ID_prob, theta_shapes, GT_prior=None):
    """
    """
    if GT_prior is None:
        GT_prior = np.ones((AD.shape[0], theta_shapes.shape[0],
                            ID_prob.shape[1])) 
        GT_prior = GT_prior / theta_shapes.shape[0]
        
    S1_gt = AD * ID_prob
    SS_gt = DP * ID_prob
    S2_gt = SS_gt - S1_gt
    
    logLik_GT = np.zeros(GT_prior.shape)
    for ig in range(logLik_GT.shape[1]):        
        logLik_GT[:, ig, :] = (S1_gt * digamma(theta_shapes[ig, 0]) +
                               S2_gt * digamma(theta_shapes[ig, 1]) - 
                               SS_gt * digamma(np.sum(theta_shapes[ig, :])))
        
    # += np.log(GT_prior)
    GT_prob = loglik_amplify(logLik_GT + np.log(GT_prior), axis=1)
    GT_prob = tensor_normalize(np.exp(GT_prob), axis=1)
    
    return GT_prob, logLik_GT



def VB_lower_bound(logLik_ID, GT_prob, ID_prob, theta_shapes, 
    theta_prior, GT_prior=None, Psi=None):
    """
    """
    if GT_prior is None:
        GT_prior = tensor_normalize(np.ones(GT_prob.shape), axis=1)
    if Psi is None:
        ID_prior = np.ones(ID_prob.shape) / ID_prob.shape[1]
    else:
        ID_prior = np.ones(ID_prob.shape) * np.log(Psi / np.sum(Psi))
        
    LB_p = np.sum(logLik_ID * ID_prob)
    LB_ID = np.sum(entropy(ID_prob.transpose(), ID_prior.transpose()))
    LB_GT = np.sum(entropy(GT_prob.transpose((1,0,2)), 
                           GT_prior.transpose((1,0,2))))
    LB_theta = beta_entropy(theta_shapes, theta_prior)
    
    # print(LB_p, LB_ID, LB_GT, LB_theta)
    return LB_p + LB_ID + LB_GT + LB_theta


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
        for ii in range(X.shape[0]):
            RV1 = (RV1 - betaln(X[ii, 0], X[ii, 1]) +
                   (X[ii, 0] - 1) * digamma(X[ii, 0]) +
                   (X[ii, 1] - 1) * digamma(X[ii, 1]) -
                   (np.sum(X[ii, :]) - 2) * digamma(np.sum(X[ii, :])))
    
    RV2 = 0
    for ii in range(X.shape[0]):
        RV2 = (RV2 - betaln(X_prior[ii, 0], X_prior[ii, 1]) +
               (X_prior[ii, 0] - 1) * digamma(X[ii, 0]) +
               (X_prior[ii, 1] - 1) * digamma(X[ii, 1]) -
               (np.sum(X_prior[ii, :]) - 2) * digamma(np.sum(X[ii, :])))
        
    return RV1 - RV2


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
        print("Randomly select 1 out %d variants" %len(idx))
        idx_use = idx[np.random.randint(len(idx))]
        
        variant_set.append(idx_use)
        barcode_set = barcode_all[idx_use]
        entropy_now = entropy_all[idx_use]
        
    if entropy_now < np.log2(K):
        print("Warning: variant_select can't distinguish all samples.")

    return entropy_now, barcode_set, variant_set
