# Core functions for Vireo model
# Author: Yuanhua Huang
# Date: 23/06/2019

# http://edwardlib.org/tutorials/probabilistic-pca
# https://github.com/allentran/pca-magic

import numpy as np
from itertools import permutations
from .vireo_base import get_ID_prob, get_GT_prob, get_theta_shapes
from .vireo_base import VB_lower_bound, tensor_normalize, loglik_amplify

def vireo_flock():
    pass

def vireo_core(AD, DP, n_donor=None, GT=None, GT_prior=None, learn_GT=True,
    theta_prior=None, learn_theta=True, Psi=None, doublet_prior=None, 
    check_doublet=True, min_iter=20, max_iter=200, epsilon_conv=1e-2,
    random_seed=None, verbose=True):
    """
    """
    if random_seed is not None:
        np.random.seed(random_seed)
        
    n_var = AD.shape[0] # n_variants and n_cells
    
    ## initialize thete
    if theta_prior is None:
        #theta_prior = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
        theta_prior = np.array([[0.1, 99.9], [50, 50], [99.9, 0.1]])
    theta_shapes = theta_prior.copy()

    ## initialize Psi
    if Psi is None:
        Psi = np.ones(n_donor) / n_donor
    else:
        Psi = Psi[:n_donor] / np.sum(Psi[:n_donor])
    ID_prob = tensor_normalize(np.random.rand(AD.shape[1], n_donor), axis=1)
    
    ## initialize GT
    if GT is None and GT_prior is None and learn_GT is False:
        print("As no GT and GT_prior is given, we change learn_GT to True.")
        learn_GT = True
    if GT_prior is None:
        GT_prior = tensor_normalize(np.ones((n_var, 3, n_donor)), axis=1)
        GT_prob, logLik_GT = get_GT_prob(AD, DP, ID_prob, 
                                         theta_shapes, GT_prior)
    else:
        # TODO: check GT_prior's shape
        GT_prob = GT_prior.copy()
    if GT is not None:
        pass
        #GT_prob = GT_to_prob(GT)

    ## VB interations
    LB = np.zeros(max_iter)
    for it in range(max_iter):
        ID_prob, GT_prob, theta_shapes, LB[it] = update_VB(AD, DP, GT_prob, 
            theta_shapes, theta_prior, GT_prior, Psi, doublet_prior,
            learn_GT=learn_GT, learn_theta=learn_theta, check_doublet=False)

        if it > min_iter:
            if LB[it] < LB[it - 1] and verbose:
                print("Warning: Lower bound decreases!\n")
            elif it == max_iter:
                print("Warning: VB did not converge!\n")
            elif LB[it] - LB[it - 1] < epsilon_conv:
                break

    ## check doublet
    if check_doublet:
        ID_prob2, GT_prob, theta_shapes, LB_doublet = update_VB(AD, DP, GT_prob, 
            theta_shapes, theta_prior, GT_prior, Psi, doublet_prior,
            learn_GT=True, learn_theta=True, check_doublet=True)

        ID_prob = ID_prob2[:, :n_donor]
        doublet_prob = ID_prob2[:, n_donor:]
    else:
        LB_doublet = None
        doublet_prob = np.zeros((ID_prob.shape[0], n_donor * (n_donor - 1) / 2))

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
        if doublet_prior is None:
            doublet_prior = 1 - GT_prob.shape[2] / GT_both.shape[2]
        Psi_both = np.ones(GT_both.shape[2]) * doublet_prior
        Psi_both[:GT_prob.shape[2]] = Psi * (1 - doublet_prior)
    else:
        Psi_both = Psi.copy()
        GT_both = GT_prob.copy()
        theta_both = theta_shapes.copy()

    ID_prob2, logLik_ID = get_ID_prob(AD, DP, GT_both, theta_both, Psi_both)
    ID_prob = ID_prob2[:, :GT_prob.shape[2]]
    
    if learn_GT:
        GT_prob, logLik_GT = get_GT_prob(AD, DP, ID_prob, 
                                         theta_shapes, GT_prior)
    if learn_theta:
        theta_shapes = get_theta_shapes(AD, DP, ID_prob, 
                                        GT_prob, theta_shapes)
    
    LB_val = VB_lower_bound(logLik_GT, GT_prob, ID_prob, theta_shapes, 
                            theta_prior, GT_prior, Psi)

    return ID_prob2, GT_prob, theta_shapes, LB_val


def add_doublet_theta(theta_shapes):
    """
    calculate theta for doublet genotype: GT=0&1 and GT=1&2 by
    averaging thire beta paramters
    
    Example
    -------
    theta_shapes = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    add_doublet_theta(theta_shapes)
    """
    #doublet GT: 0_1 and 1_2
    theta_shapes2 = np.zeros((2,2)) 
    for ii in range(2):
        theta_use  = theta_shapes[ii:(ii + 2), :]
        theta_mean = np.mean(tensor_normalize(theta_use, axis=1), axis=0)
        shape_sum  = np.sqrt(np.sum(theta_use[0, :]) * np.sum(theta_use[1, :]))
        theta_shapes2[ii, :] = theta_mean * shape_sum
    return np.append(theta_shapes, theta_shapes2, axis=0)


def add_doublet_GT(GT_prob):
    """
    Add doublet genotype by summarizing their probability:
    New GT has five categories: 0, 1, 2, 1.5, 2.5
    """
    db_idx = np.array(list(permutations(range(GT_prob.shape[2]), 2)))
    idx1 = db_idx[:, 0]
    idx2 = db_idx[:, 1]
    
    ## GT_prob has three genotypes: 0, 1, 2;
    GT_prob2 = np.zeros((GT_prob.shape[0], 5, db_idx.shape[0]))
    
    GT_prob2[:, 0, :] = (GT_prob[:, 0, idx1] * GT_prob[:, 0, idx2])
    GT_prob2[:, 1, :] = (GT_prob[:, 0, idx1] * GT_prob[:, 2, idx2] +
                         GT_prob[:, 1, idx1] * GT_prob[:, 1, idx2] +
                         GT_prob[:, 2, idx1] * GT_prob[:, 0, idx2])
    GT_prob2[:, 2, :] = (GT_prob[:, 2, idx1] * GT_prob[:, 2, idx2])
    GT_prob2[:, 3, :] = (GT_prob[:, 0, idx1] * GT_prob[:, 1, idx2] + 
                         GT_prob[:, 1, idx1] * GT_prob[:, 0, idx2])
    GT_prob2[:, 4, :] = (GT_prob[:, 1, idx1] * GT_prob[:, 2, idx2] + 
                         GT_prob[:, 2, idx1] * GT_prob[:, 1, idx2])
    
    GT_prob2 = tensor_normalize(GT_prob2, axis=1)
    GT_prob1 = np.append(GT_prob, np.zeros((GT_prob.shape[0], 2, 
                                            GT_prob.shape[2])), axis=1)
    return np.append(GT_prob1, GT_prob2, axis=2)
    
