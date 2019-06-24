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
    
    logLik_ID += np.log(Psi / np.sum(Psi))
    logLik_vec = logsumexp(logLik_ID, axis=1)
    
    ID_prob = np.exp(loglik_amplify(logLik_ID, axis=1))
    ID_prob = tensor_normalize(ID_prob, axis=1)
    
    return ID_prob, np.sum(logLik_vec)
    

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



def VB_lower_bound(logLik_GT, GT_prob, ID_prob, theta_shapes, 
    theta_prior, GT_prior=None, Psi=None):
    """
    """
    if GT_prior is None:
        GT_prior = tensor_normalize(np.ones(GT_prob.shape), axis=1)
    if Psi is None:
        ID_prior = np.ones(ID_prob.shape) / ID_prob.shape[1]
    else:
        ID_prior = np.ones(ID_prob.shape) * np.log(Psi / np.sum(Psi))
        
    LB_p = np.sum(logLik_GT * GT_prob) #+ sum(W_vec)
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
