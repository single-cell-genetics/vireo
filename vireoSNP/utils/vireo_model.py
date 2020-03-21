# Core functions for Vireo model
# Author: Yuanhua Huang
# Date: 30/08/2019

# http://edwardlib.org/tutorials/probabilistic-pca
# https://github.com/allentran/pca-magic

import sys
import itertools
import numpy as np
from scipy.stats import entropy
from scipy.special import digamma
from .vireo_base import normalize, loglik_amplify, beta_entropy

def vireo_core(AD, DP, n_donor=None, GT_prior=None, learn_GT=True,
    theta_prior=None, learn_theta=True, ASE_mode=False, 
    Psi=None, ID_prob_init=None, doublet_prior=None, check_doublet=True, 
    min_iter=20, max_iter=100, min_GP=0.00001, epsilon_conv=1e-2, 
    random_seed=None, verbose=False):
    """
    Vireo core function to cluster the cells into donors.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    if n_donor is None:
        if len(GT_prior.shape) < 3 or GT_prior.shape[1] < 2:
            print("Error: no n_donor and GT_prior has < 2 donors.")
            sys.exit(1)
        else:
            n_donor = GT_prior.shape[1]
        
    n_var = AD.shape[0] # n_variants
    
    ## initialize thete
    if theta_prior is None:
        #theta_prior = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
        theta_prior = np.array([[0.1, 99.9], [50, 50], [99.9, 0.1]])
    theta_shapes = theta_prior.copy()
    if ASE_mode and len(theta_prior.shape) == 2:
        theta_prior = np.repeat(np.expand_dims(theta_prior, 2), n_var, axis=2)
        theta_shapes = np.repeat(np.expand_dims(theta_shapes, 2), n_var, axis=2)
    n_gt = theta_shapes.shape[0] # number of genotype categories

    ## initialize Psi
    if Psi is None:
        Psi = np.ones(n_donor) / n_donor
    else:
        Psi = Psi[:n_donor] / np.sum(Psi[:n_donor])
    if ID_prob_init is None:
        ID_prob = normalize(np.random.rand(AD.shape[1], n_donor))
    else:
        ID_prob = normalize(ID_prob_init.copy())
    
    ## initialize GT
    if GT_prior is None:
        GT_prior = normalize(np.ones((n_var, n_donor, n_gt)))
        GT_prob, logLik_GT = get_GT_prob(AD, DP, ID_prob, 
                                         theta_shapes, GT_prior)
        if learn_GT is False:
            print("As GT_prior is not given, we change learn_GT to True.")
            learn_GT = True
    else:
        GT_prob = GT_prior.copy()
        GT_prior[GT_prior < min_GP] = min_GP
        GT_prior[GT_prior > 1 - min_GP] = 1 - min_GP
        GT_prior = normalize(GT_prior)

    #TODO: check if there is a better way to deal with GT imcompleteness
    if GT_prior.shape[1] < n_donor:
        _add_n = n_donor - GT_prior.shape[1]
        GT_prior = np.append(GT_prior, 
            normalize(np.ones((n_var, n_gt, _add_n)), axis=1))
        GT_prob = GT_prior.copy()
        if learn_GT is False:
            print("As GT_prior is not complete, we change learn_GT to True.")
            learn_GT = True
    elif GT_prior.shape[1] > n_donor:
        print("Warning: n_donor is smaller than samples in GT_prior, hence we "
              "ignore n_donor.")
        n_donor = GT_prior.shape[1]

    # check if n_gt is matched to GT_prior
    if GT_prior.shape[2] != n_gt:
        print("Error: number of GT categories not matched: theta and GT_prior")
        sys.exit(1)

    ## VB interations
    LB = np.zeros(max_iter)
    for it in range(max_iter):
        ID_prob, GT_prob, theta_shapes, LB[it] = update_VB(AD, DP, GT_prob, 
            theta_shapes, theta_prior, GT_prior, Psi, doublet_prior,
            learn_GT=learn_GT, learn_theta=learn_theta, 
            check_doublet=check_doublet)

        if it > min_iter:
            if LB[it] < LB[it - 1]:
                if verbose:
                    print("Warning: Lower bound decreases!\n")
            elif it == max_iter - 1: 
                if verbose:
                    print("Warning: VB did not converge!\n")
            elif LB[it] - LB[it - 1] < epsilon_conv:
                break

    ## one-off check doublet
    if check_doublet:
        ID_prob2, GT_prob, theta_shapes, LB_doublet = update_VB(AD, DP, GT_prob, 
            theta_shapes, theta_prior, GT_prior, Psi, doublet_prior,
            learn_GT=True, learn_theta=learn_theta, check_doublet=True)

        ID_prob = ID_prob2[:, :n_donor]
        doublet_prob = ID_prob2[:, n_donor:]
    else:
        LB_doublet = LB[it]
        n_donor_doublt = int(n_donor * (n_donor - 1) / 2)
        doublet_prob = np.zeros((ID_prob.shape[0], n_donor_doublt))

    RV = {}
    RV['ID_prob'] = ID_prob
    RV['GT_prob'] = GT_prob
    RV['doublet_prob'] = doublet_prob
    RV['theta_shapes'] = theta_shapes
    RV['LB_list'] = LB[: it+1]
    RV['LB_doublet'] = LB_doublet
    return RV


def update_VB(AD, DP, GT_prob, theta_shapes, theta_prior, GT_prior, 
    Psi, doublet_prior=None, learn_GT=True, learn_theta=True, 
    check_doublet=False):
    """
    Update the parameters of each component of the variantional 
    distribution. 
    
    The doublet probability can be created by doublet genotypes
    """
    if check_doublet:
        GT_both = add_doublet_GT(GT_prob)
        theta_both = add_doublet_theta(theta_shapes)
        n_doublet_pair = GT_both.shape[1] - GT_prob.shape[1]
        if doublet_prior is None:
            doublet_prior = min(0.5, AD.shape[1] / 100000)
            
        Psi_both = np.append(Psi * (1 - doublet_prior), 
                             (np.ones(n_doublet_pair) / n_doublet_pair * 
                              doublet_prior))
    else:
        Psi_both = Psi.copy()
        GT_both = GT_prob.copy()
        theta_both = theta_shapes.copy()

    ID_prob2, logLik_ID = get_ID_prob(AD, DP, GT_both, theta_both, Psi_both)
    ID_prob = ID_prob2[:, :GT_prob.shape[1]]
    
    if learn_GT:
        GT_prob, logLik_GT = get_GT_prob(AD, DP, ID_prob, 
                                         theta_shapes, GT_prior)
    if learn_theta:
        theta_shapes = get_theta_shapes(AD, DP, ID_prob, 
                                        GT_prob, theta_prior)
    
    ### check how to calculate lower bound for when detecting doublets
    LB_val = VB_lower_bound(logLik_ID, GT_prob, ID_prob2, theta_shapes, 
                            theta_prior, GT_prior, Psi_both)

    return ID_prob2, GT_prob, theta_shapes, LB_val



def get_theta_shapes(AD, DP, ID_prob, GT_prob, theta_prior):
    """
    """
    S1_gt = AD * ID_prob
    SS_gt = DP * ID_prob
    S2_gt = SS_gt - S1_gt
    
    theta_shapes = theta_prior.copy()
    for ig in range(theta_shapes.shape[0]):
        _axis = 1 if len(theta_shapes.shape) == 3 else None
        theta_shapes[ig, 0] += np.sum(S1_gt * GT_prob[:, :, ig], axis=_axis)
        theta_shapes[ig, 1] += np.sum(S2_gt * GT_prob[:, :, ig], axis=_axis)
    return theta_shapes

def get_ID_prob(AD, DP, GT_prob, theta_shapes, Psi=None):
    """
    """
    if Psi is None:
        Psi = np.ones(GT_prob.shape[1]) / GT_prob.shape[1]

    BD = DP - AD
    logLik_ID = np.zeros((AD.shape[1], GT_prob.shape[1]))
    for ig in range(GT_prob.shape[2]):
        _digmma1 = digamma(theta_shapes[ig, 0]).reshape(-1, 1)
        _digmma2 = digamma(theta_shapes[ig, 1]).reshape(-1, 1)
        _digmmas = digamma(theta_shapes[ig, :].sum(axis=0)).reshape(-1, 1)
        S1 = AD.transpose() * (GT_prob[:, :, ig] * _digmma1)
        S2 = BD.transpose() * (GT_prob[:, :, ig] * _digmma2)
        SS = DP.transpose() * (GT_prob[:, :, ig] * _digmmas)
        logLik_ID += (S1 + S2 - SS)
    
    Psi_norm = np.log(Psi / np.sum(Psi))
    ID_prob = np.exp(loglik_amplify(logLik_ID + Psi_norm, axis=1))
    ID_prob = normalize(ID_prob, axis=1)
    
    return ID_prob, logLik_ID
    

def get_GT_prob(AD, DP, ID_prob, theta_shapes, GT_prior=None):
    """
    """
    if GT_prior is None:
        GT_prior = np.ones((AD.shape[0], ID_prob.shape[1], 
                            theta_shapes.shape[0])) 
        GT_prior = GT_prior / theta_shapes.shape[0]
        
    S1_gt = AD * ID_prob
    SS_gt = DP * ID_prob
    S2_gt = SS_gt - S1_gt
    
    logLik_GT = np.zeros(GT_prior.shape)
    for ig in range(logLik_GT.shape[2]):        
        _digmma1 = digamma(theta_shapes[ig, 0]).reshape(-1, 1)
        _digmma2 = digamma(theta_shapes[ig, 1]).reshape(-1, 1)
        _digmmas = digamma(theta_shapes[ig, :].sum(axis=0)).reshape(-1, 1)
        logLik_GT[:, :, ig] = (S1_gt * _digmma1 + 
                               S2_gt * _digmma2 - 
                               SS_gt * _digmmas)
        
    # += np.log(GT_prior)
    GT_prob = loglik_amplify(logLik_GT + np.log(GT_prior), axis=2)
    GT_prob = normalize(np.exp(GT_prob), axis=2)
    
    return GT_prob, logLik_GT



def VB_lower_bound(logLik_ID, GT_prob, ID_prob, theta_shapes, 
    theta_prior, GT_prior=None, Psi=None):
    """
    """
    if GT_prior is None:
        GT_prior = normalize(np.ones(GT_prob.shape), axis=2)
    if Psi is None:
        ID_prior = np.ones(ID_prob.shape) / ID_prob.shape[1]
    else:
        ID_prior = np.ones(ID_prob.shape) * np.log(Psi / np.sum(Psi))
        
    LB_p = np.sum(logLik_ID * ID_prob)
    KL_ID = -np.sum(entropy(ID_prob, ID_prior, axis=1))
    KL_GT = -np.sum(entropy(GT_prob, GT_prior, axis=2))
    KL_theta = -beta_entropy(theta_shapes, theta_prior)
    
    # print(LB_p, KL_ID, KL_GT, KL_theta)
    return LB_p - KL_ID - KL_GT - KL_theta


def add_doublet_theta(theta_shapes):
    """
    calculate theta for doublet genotype: GT=0&1, GT=0&2, and GT=1&2 by
    averaging thire beta paramters
    
    Example
    -------
    theta_shapes = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    add_doublet_theta(theta_shapes)
    """
    # TODO: support reduced GT for relatives
    combn_iter = itertools.combinations(range(theta_shapes.shape[0]), 2)
    db_idx = np.array([x for x in combn_iter])

    _theta_p1 = theta_shapes[db_idx[:, 0]]
    _theta_p2 = theta_shapes[db_idx[:, 1]]

    _theta_mean = (normalize(_theta_p1, axis=1) + 
                   normalize(_theta_p2, axis=1)) / 2.0
    _theta_sum  = np.sqrt(np.sum(_theta_p1, axis=1, keepdims=True) * 
                          np.sum(_theta_p2, axis=1, keepdims=True))
    
    theta_shapes_db = _theta_mean * _theta_sum

    return np.append(theta_shapes, theta_shapes_db, axis=0)


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
