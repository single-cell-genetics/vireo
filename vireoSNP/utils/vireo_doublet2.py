## Prediction doublets

import itertools
import numpy as np
import multiprocessing
from scipy.stats import entropy
from scipy.sparse import csc_matrix
from scipy.special import logsumexp, digamma, betaln
from .vireo_base import normalize, loglik_amplify
from .variant_select import variant_ELBO_gain


def predict_doublet2(vobj, AD, DP, update_GT=True, update_ID=True, 
    doublet_rate_prior=0.5, doublet_resolution=0.1):
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
    def _get_theta_size(ID_prob, AD, DP, cell_use):
        """Coordinate ascent for updating theta posterior parameters
        """
        BD = DP - AD
        _theta_s1 = AD[:, cell_use] @ ID_prob[cell_use, :]  #(n_var, n_donor)
        _theta_s2 = BD[:, cell_use] @ ID_prob[cell_use, :]  #(n_var, n_donor)
        # _theta_s1 += vobj.theta_s1_prior
        # _theta_s2 += vobj.theta_s2_prior

        # _theta_s1 += 0.5
        # _theta_s2 += 0.5

        _beta_mu = _theta_s1 / (_theta_s1 + _theta_s2)
        _beta_sum = _theta_s1 + _theta_s2

        return _beta_mu, _beta_sum

    def _get_ID_prob(theta_s1_both, theta_s2_both, AD, DP, snp_use, 
        ID_prior_both=None,
        mode='binom'):
        """Calculate cell assignment probability
        """
        BD = DP - AD

        AD = AD[snp_use, :]
        BD = BD[snp_use, :]
        DP = DP[snp_use, :]
        theta_s1_both = theta_s1_both[snp_use, :]
        theta_s2_both = theta_s2_both[snp_use, :]
    
        if mode == 'betabin':
            _E_logLik_mat = np.zeros((AD.shape[1], theta_s2_both.shape[1]))
            for i in range(AD.shape[1]):
                _E_logLik_mat[i, :] = (
                    betaln(AD[:, i].A + theta_s1_both, BD[:, i].A + theta_s2_both) - 
                    betaln(theta_s1_both, theta_s2_both)
                ).sum(axis=0)
        else:
            _E_logLik_mat = (
                AD.T @ digamma(theta_s1_both) + 
                BD.T @ digamma(theta_s2_both) -
                DP.T @ digamma(theta_s1_both + theta_s2_both)
            ) # shape: (n_cell, n_combined_donor)

        if ID_prior_both is None:
            ID_prior_both = np.ones(_E_logLik_mat.shape)

        _ID_prob_both = normalize(np.exp(loglik_amplify(
            _E_logLik_mat + np.log(ID_prior_both))))

        return _ID_prob_both

    doublet_prob_sum = np.zeros(AD.shape[1])
    ID_prob = vobj.ID_prob.copy()

    _db_props = np.arange(
        doublet_resolution, (1 - doublet_resolution + 0.01),
        doublet_resolution)

    # _db_props = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    _db_props = [0.3, 0.4, 0.5, 0.6, 0.7]

    for _iter in range(5):
        _cell_idx = doublet_prob_sum < 0.9
        _beta_mu, _beta_sum = _get_theta_size(ID_prob, AD, DP, _cell_idx)

        # _snp_idx = ((np.median(_beta_sum, axis=1) >= 5) * 
        #             (np.var(_beta_mu, axis=1) >= 0.03))
        
        min_ELBO_gain = AD.shape[1] / 100.0
        _ELBO_gain = variant_ELBO_gain(ID_prob, AD, DP)
        _snp_idx = _ELBO_gain >= np.sqrt(AD.shape[1]) / 3.0
        print("%d out %d SNPs selected for doublet detection: ELBO_gain > %.1f" 
              %(sum(_snp_idx), len(_snp_idx), min_ELBO_gain))

        beta_mu_both, beta_sum_both = add_doublet_theta(
            _beta_mu, _beta_sum, _db_props)

        theta_s1_both = beta_mu_both * beta_sum_both
        theta_s2_both = beta_sum_both - theta_s1_both

        # theta_s1_both[:, beta_mu_both.shape[1]:] += 1
        # theta_s2_both[:, beta_mu_both.shape[1]:] += 1
        theta_s1_both += 0.1
        theta_s2_both += 0.1

        n_doublet_pair = beta_mu_both.shape[1] - vobj.n_donor
        if doublet_rate_prior is None:
            doublet_rate_prior = min(0.5, AD.shape[1] / 100000)
            
        ID_prior_both = np.append(
            vobj.ID_prior * (1 - doublet_rate_prior), 
            np.ones((vobj.n_cell, n_doublet_pair)) / n_doublet_pair * 
            doublet_rate_prior, axis=1)

        # Calculate assignment probability (same as update_ID_prob())
        ID_prob_both = _get_ID_prob(theta_s1_both, theta_s2_both, AD, DP, 
                                    snp_use=_snp_idx, 
                                    ID_prior_both=ID_prior_both, mode='binom')

        _doublet_prob_sum = ID_prob_both[:, vobj.n_donor:].sum(1)
        # print(np.sum((doublet_prob_sum > 0.5)), np.sum((_doublet_prob_sum > 0.5)))

        if np.mean((_doublet_prob_sum > 0.5) == (doublet_prob_sum > 0.5)) == 1:
            pass
            # break
        else:
            ID_prob = normalize(ID_prob_both[:, :vobj.n_donor])
            doublet_prob_sum = _doublet_prob_sum + 0.0

    doublet_prob = np.max(ID_prob_both[:, vobj.n_donor:].reshape(
        ID_prob_both.shape[0], len(_db_props), -1), 1)

    return doublet_prob, ID_prob_both[:, :vobj.n_donor]


def add_doublet_theta(beta_mu, beta_sum, proportions=[0.5]):
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

    # beta_mu_db = (beta_mu[:, db_idx[:, 0]] + beta_mu[:, db_idx[:, 1]]) / 2.0
    # # beta_sum_db = np.sqrt(beta_sum[:, db_idx[:, 0]] * beta_sum[:, db_idx[:, 1]])
    # beta_sum_db = np.dstack(
    #     (beta_sum[:, db_idx[:, 0]], beta_sum[:, db_idx[:, 1]])).mean(2)

    beta_s1 = beta_mu * beta_sum
    beta_s2 = (1 - beta_mu) * beta_sum
    beta_s1_db = np.hstack(
        [x * beta_s1[:, db_idx[:, 0]] + (1 - x) * beta_s1[:, db_idx[:, 1]] \
         for x in proportions]
    )
    beta_s2_db = np.hstack(
        [x * beta_s2[:, db_idx[:, 0]] + (1 - x) * beta_s2[:, db_idx[:, 1]] \
         for x in proportions]
    )
    
    # beta_s1_db = np.sqrt(beta_s1[:, db_idx[:, 0]] * beta_s1[:, db_idx[:, 1]])
    # beta_s2_db = np.sqrt(beta_s2[:, db_idx[:, 0]] * beta_s2[:, db_idx[:, 1]])
    
    beta_mu_db = beta_s1_db / (beta_s1_db + beta_s2_db)
    beta_sum_db = beta_s1_db + beta_s2_db

    return (np.append(beta_mu, beta_mu_db, axis=-1), 
            np.append(beta_sum, beta_sum_db, axis=-1))

