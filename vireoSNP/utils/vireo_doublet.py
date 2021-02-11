## Prediction doublets

import itertools
import numpy as np
import multiprocessing
from scipy.stats import entropy
from scipy.sparse import csc_matrix
from scipy.special import logsumexp, digamma, betaln
from .vireo_base import normalize, loglik_amplify

def predict_doublet(vobj, AD, DP, update_GT=True, update_ID=True, 
    doublet_rate_prior=None):
    """Predict doublet with fitted Vireo model

    Parameters
    ----------
    vobj : Vireo object
        Fitted Vireo object before predicting doublets
    AD : scipy.sparse.csc_matrix (n_var, n_cell)
        Sparse count matrix for alternative allele
    DP : scipy.sparse.csc_matrix (n_var, n_cell)
        Sparse count matrix for depths, alternative + refeerence alleles
    update_GT : bool
        Whether updating GT_prob after removing doublet_prob
    update_GT : bool
        Whether updating ID_prob by removing doublet_prob
    doublet_rate_prior : float
        Prior value of doublet rate

    Returns
    -------
    A tuple of two numpy arrays (doublet_prob, ID_prob)

    doublet_prob : numpy array (n_cell, n_donor * (n_donor - 1) / 2)
        Assignment probability of a cells to any doublet (donor pair)
    ID_prob : numpy array (n_cell, n_donor)
        updated ID_prob by removing doublet_prob
    """
    GT_both = add_doublet_GT(vobj.GT_prob)
    beta_mu_both, beta_sum_both = add_doublet_theta(vobj.beta_mu, 
                                                    vobj.beta_sum)

    n_doublet_pair = GT_both.shape[1] - vobj.GT_prob.shape[1]
    if doublet_rate_prior is None:
        doublet_rate_prior = min(0.5, AD.shape[1] / 100000)
        
    ID_prior_both = np.append(
        vobj.ID_prior * (1 - doublet_rate_prior), 
        np.ones((vobj.n_cell, n_doublet_pair)) / n_doublet_pair * 
        doublet_rate_prior, axis=1)

    # Calculate assignment probability (same as update_ID_prob())
    BD = DP - AD
    logLik_ID = np.zeros((AD.shape[1], GT_both.shape[1]))
    _digamma1 = np.expand_dims(digamma(beta_sum_both * beta_mu_both), 1)
    _digamma2 = np.expand_dims(digamma(beta_sum_both * (1 - beta_mu_both)), 1)
    _digammas = np.expand_dims(digamma(beta_sum_both), 1)
    for ig in range(GT_both.shape[2]):
        S1 = AD.T @ (GT_both[:, :, ig] * _digamma1[:, :, ig])
        S2 = BD.T @ (GT_both[:, :, ig] * _digamma2[:, :, ig])
        SS = DP.T @ (GT_both[:, :, ig] * _digammas[:, :, ig])
        logLik_ID += (S1 + S2 - SS)

    ID_prob_both = normalize(np.exp(loglik_amplify(
        logLik_ID + np.log(ID_prior_both))))

    if update_ID:
        vobj.ID_prob = ID_prob_both[:, :vobj.n_donor]

    if update_GT:
        if update_ID:
            vobj.update_GT_prob(AD, DP)
        else:
            print("For update_GT, please turn on update_ID.")
    
    return ID_prob_both[:, vobj.n_donor:], ID_prob_both[:, :vobj.n_donor]


def add_doublet_theta(beta_mu, beta_sum):
    """
    calculate theta for doublet genotype: GT=0&1, GT=0&2, and GT=1&2 by
    averaging thire beta paramters
    
    Example
    -------
    add_doublet_theta(np.array([[0.01, 0.5, 0.99]]), np.array([[30, 6, 30]]))
    """
    # TODO: support reduced GT for relatives
    combn_iter = itertools.combinations(range(beta_mu.shape[1]), 2)
    db_idx = np.array([x for x in combn_iter])

    beta_mu_db = (beta_mu[:, db_idx[:, 0]] + beta_mu[:, db_idx[:, 1]]) / 2.0
    beta_sum_db = np.sqrt(beta_sum[:, db_idx[:, 0]] * beta_sum[:, db_idx[:, 1]])

    return (np.append(beta_mu, beta_mu_db, axis=-1), 
            np.append(beta_sum, beta_sum_db, axis=-1))


def add_doublet_GT(GT_prob):
    """
    Add doublet genotype by summarizing their probability:
    New GT has five categories: 0, 1, 2, 1.5, 2.5
    TODO: New GT has six categories: 0, 1, 2, 0_1, 0_2, 1_2
    """
    combn_iter = itertools.combinations(range(GT_prob.shape[2]), 2)
    gt_idx = np.array([x for x in combn_iter]) # GT combination
    g_idx1 = gt_idx[:, 0]
    g_idx2 = gt_idx[:, 1]

    combn_iter = itertools.combinations(range(GT_prob.shape[1]), 2)
    sp_idx = np.array([x for x in combn_iter]) # sample combination
    s_idx1 = sp_idx[:, 0]
    s_idx2 = sp_idx[:, 1]
    
    ## GT_prob has three genotypes: 0, 1, 2;
    n_gt = GT_prob.shape[2]
    GT_prob2 = np.zeros((GT_prob.shape[0], sp_idx.shape[0],
                         n_gt + gt_idx.shape[0]))

    GT_prob2[:, :, :n_gt] = (GT_prob[:, s_idx1, :] * 
                             GT_prob[:, s_idx2, :])
    GT_prob2[:, :, n_gt:] = (GT_prob[:, s_idx1, :][:, :, g_idx1] * 
                             GT_prob[:, s_idx2, :][:, :, g_idx2] +
                             GT_prob[:, s_idx1, :][:, :, g_idx2] * 
                             GT_prob[:, s_idx2, :][:, :, g_idx1])
    
    GT_prob2 = normalize(GT_prob2, axis=2)
    GT_prob1 = np.append(GT_prob, 
        np.zeros((GT_prob.shape[0], GT_prob.shape[1], gt_idx.shape[0])), axis=2)
    return np.append(GT_prob1, GT_prob2, axis=1)


def _fit_EM_ambient(AD, DP, theta_mat, n_donor, 
    max_iter=200, min_iter=5, epsilon_conv=1e-3, verbose=False):
    """Estimate ambient RNA abundance by EM algorithm
    """
    BD = DP - AD
    psi = np.random.dirichlet([1] * n_donor)
    logLik = np.zeros(max_iter)
    for it in range(max_iter):
        # E step: expectation of count assignment probability
        Z1 = theta_mat * np.expand_dims(psi, 0)
        Z1 = Z1 / np.sum(Z1, axis=1, keepdims=True)
        
        Z0 = (1 - theta_mat) * np.expand_dims(psi, 0)
        Z0 = Z0 / np.sum(Z0, axis=1, keepdims=True)

        # M step: maximize logLikehood over psi and theta
        psi_raw = np.dot(AD, Z1) + np.dot(BD, Z0)
        psi = psi_raw / np.sum(psi_raw)
        
        # Likelihood and check convergence
        theta_vct = np.dot(theta_mat, psi)
        logLik[it] = np.sum(
            AD * np.log(theta_vct) + BD * np.log(1 - theta_vct))
        if it > min_iter:
            if logLik[it] < logLik[it - 1]:
                if verbose:
                    print("Warning: logLikelihood decreases!\n")
            elif it == max_iter - 1:
                if verbose:
                    print("Warning: VB did not converge!\n")
            elif logLik[it] - logLik[it - 1] < epsilon_conv:
                break
    logLik_RV = logLik[:it]
    return psi
    

def predit_ambient(vobj, AD, DP, nproc=10):
    """Predict fraction of ambient RNA contaimination.
    Still under development
    """
    ### detect ambient RNA for each cell
    import timeit
    start = timeit.default_timer()

    # theta_mat = (AD @ vobj.ID_prob + 0.1) / (DP @ vobj.ID_prob + 0.2)
    theta_mat = np.tensordot(vobj.GT_prob, vobj.beta_mu[0, :], axes=(2, 0))

    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes = nproc)
        for i in range(AD.shape[1]): 
            _ad = AD[:, i].toarray().reshape(-1)
            _dp = DP[:, i].toarray().reshape(-1)
            result.append(pool.apply_async(_fit_EM_ambient, 
            (_ad, _dp, theta_mat, vobj.n_donor), callback = None))
        pool.close()
        pool.join()
        Psi_mat = np.array([res.get() for res in result])
    else:
        Psi_mat = np.zeros((AD.shape[1], vobj.n_donor))
        for i in range(AD.shape[1]):
            _ad = AD[:, i].toarray().reshape(-1)
            _dp = DP[:, i].toarray().reshape(-1)
            Psi_mat[i, :] = _fit_EM_ambient(
                _ad, _dp, theta_mat, vobj.n_donor)

        return Psi_mat

    stop = timeit.default_timer()
    print('Time: ', stop - start)

    return Psi_mat

