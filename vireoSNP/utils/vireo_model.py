# Core functions for Vireo model
# Author: Yuanhua Huang
# Date: 30/08/2019

# http://edwardlib.org/tutorials/probabilistic-pca
# https://github.com/allentran/pca-magic

import sys
import itertools
import numpy as np
import multiprocessing
from .vireo_base import get_ID_prob, get_GT_prob, get_theta_shapes
from .vireo_base import VB_lower_bound, tensor_normalize, loglik_amplify

def show_progress(RV=None):
    return RV

def vireo_flock(AD, DP, GT_prior=None, n_donor=None, n_extra_donor=2, 
                extra_donor_mode="distance", n_init=20, random_seed=None, 
                check_doublet=True, **kwargs):
    """
    A wrap function to run vireo twice, with the first step searching many 
    initialations.
    """
    ## random seed setting
    if random_seed is not None:
        np.random.seed(random_seed)

    if n_donor is None:
        if GT_prior is None:
            print("[vireo] Error: requiring n_donor or GT_prior.")
            sys.exit()
        else:
            n_donor = GT_prior.shape[2]

    ## warm initialization optionally with extra components
    _ID_prob = None
    if GT_prior is None or n_donor > GT_prior.shape[2]:
        n_donor_run1 = int(n_donor + n_extra_donor)
        print("[vireo] warm-up: %d random initializations for %d clusters..." 
              %(n_init, n_donor_run1))

        ID_prob_list = []
        for i in range(n_init):
            _ID_prob = np.random.rand(AD.shape[1], n_donor_run1)
            ID_prob_list.append(tensor_normalize(_ID_prob, axis=1))

        result = []
        for _ID_prob in ID_prob_list:
            result.append(vireo_core(AD, DP, n_donor=n_donor_run1,
                GT_prior = None, ID_prob_init=_ID_prob, min_iter=5, max_iter=15, 
                verbose=False, check_doublet=False, **kwargs))

        LB_list = [x['LB_list'][-1] for x in result]
        res1 = result[np.argmax(LB_list)]

        _ID_prob = donor_select(res1['GT_prob'], res1['ID_prob'], n_donor, 
                                mode=extra_donor_mode)

        print("[vireo] warm-up: lower bound ranges [%.1f, %.1f, %.1f]" 
              %(min(LB_list), np.median(LB_list), max(LB_list)))
    else:
        _ID_prob = None

    ## pre-run: tune the genotype prior (if genotype prior is uncertain)
    GT_prior_use = GT_prior
    if GT_prior is not None and n_donor < GT_prior.shape[2]:
        print("[vireo] pre-RUN: finding %d from %d donors with GT ..." 
              %(n_donor, GT_prior.shape[2]))
        res1 = vireo_core(AD, DP, GT_prior = GT_prior, n_donor=None, 
                          ID_prob_init=_ID_prob, check_doublet=False, **kwargs)
        _donor_cnt = np.sum(res1['ID_prob'], axis=0)
        _donor_idx = np.argsort(_donor_cnt)[::-1]
        GT_prior_use = GT_prior[:, :, _donor_idx[:n_donor]]

        print("\t".join(["donor%d" %x for x in _donor_idx]))
        print("\t".join(["%.0f" %_donor_cnt[x] for x in _donor_idx]))

    elif GT_prior is not None and n_donor > GT_prior.shape[2]:
        print("[vireo] pre-RUN: finding %d from %d donors without GT ..." 
              %(n_donor - GT_prior.shape[2], n_donor))
        res1 = vireo_core(AD, DP, GT_prior=None, n_donor=n_donor, 
                          ID_prob_init=_ID_prob, check_doublet=check_doublet, 
                          **kwargs)
        GT_prior_use = res1['GT_prob']
        idx = greed_match(GT_prior, GT_prior_use)
        GT_prior_use[:, :, idx] = GT_prior

        _donor_cnt = np.sum(res1['ID_prob'], axis=0) 
        idx_ordered = np.append(idx, np.delete(np.arange(n_donor), idx))
        print("\t".join(["donor%d" %x for x in idx_ordered]))
        print("\t".join(["%.0f" %_donor_cnt[x] for x in idx_ordered]))

    ## main run
    print("[vireo] main RUN with warm initials and tuned GT ...")
    res1 = vireo_core(AD, DP, GT_prior=GT_prior_use, n_donor=n_donor, 
                      ID_prob_init=_ID_prob, check_doublet=check_doublet, 
                      **kwargs)
    print("[vireo] main RUN: %d iterations; lower bound %.1f" 
          %(len(res1['LB_list']), res1['LB_list'][-1]))

    ## print the beta parameters
    print("[vireo] beta parameters for binomial rate:")
    np.set_printoptions(formatter={'float': lambda x: format(x, '.1f')})
    print(res1['theta_shapes'])

    return res1


def vireo_core(AD, DP, n_donor=None, GT_prior=None, learn_GT=True,
    theta_prior=None, learn_theta=True, Psi=None, ID_prob_init=None, 
    doublet_prior=None, check_doublet=True, min_iter=20, max_iter=100, 
    min_GP=0.00001, epsilon_conv=1e-2, random_seed=None, verbose=False):
    """
    Vireo core function to cluster the cells into donors.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    if n_donor is None:
        if len(GT_prior.shape) < 3 or GT_prior.shape[2] < 2:
            print("Error: no n_donor and GT_prior has < 2 donors.")
            sys.exit(1)
        else:
            n_donor = GT_prior.shape[2]
        
    n_var = AD.shape[0] # n_variants
    
    ## initialize thete
    if theta_prior is None:
        #theta_prior = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
        theta_prior = np.array([[0.1, 99.9], [50, 50], [99.9, 0.1]])
    theta_shapes = theta_prior.copy()
    n_gt = theta_shapes.shape[0] # number of genotype categories

    ## initialize Psi
    if Psi is None:
        Psi = np.ones(n_donor) / n_donor
    else:
        Psi = Psi[:n_donor] / np.sum(Psi[:n_donor])
    if ID_prob_init is None:
        ID_prob = tensor_normalize(np.random.rand(AD.shape[1], n_donor), axis=1)
    else:
        ID_prob = tensor_normalize(ID_prob_init.copy(), axis=1)
    
    ## initialize GT
    if GT_prior is None:
        GT_prior = tensor_normalize(np.ones((n_var, n_gt, n_donor)), axis=1)
        GT_prob, logLik_GT = get_GT_prob(AD, DP, ID_prob, 
                                         theta_shapes, GT_prior)
        if learn_GT is False:
            print("As GT_prior is not given, we change learn_GT to True.")
            learn_GT = True
    else:
        GT_prob = GT_prior.copy()
        GT_prior[GT_prior < min_GP] = min_GP
        GT_prior[GT_prior > 1 - min_GP] = 1 - min_GP
        GT_prior = tensor_normalize(GT_prior, axis=1)

    #TODO: check if there is a better way to deal with GT imcompleteness
    if GT_prior.shape[2] < n_donor:
        _add_n = n_donor - GT_prior.shape[2]
        GT_prior = np.append(GT_prior, 
            tensor_normalize(np.ones((n_var, n_gt, _add_n)), axis=1), axis=2)
        GT_prob = GT_prior.copy()
        if learn_GT is False:
            print("As GT_prior is not complete, we change learn_GT to True.")
            learn_GT = True
    elif GT_prior.shape[2] > n_donor:
        print("Warning: n_donor is smaller than samples in GT_prior, hence we "
              "ignore n_donor.")
        n_donor = GT_prior.shape[2]

    # check if n_gt is matched to GT_prior
    if GT_prior.shape[1] != n_gt:
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
        n_doublet_pair = GT_both.shape[2] - GT_prob.shape[2]
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
    ID_prob = ID_prob2[:, :GT_prob.shape[2]]
    
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

    theta_shapes2 = np.zeros((db_idx.shape))
    for ii in range(theta_shapes2.shape[0]):
        theta_use  = theta_shapes[db_idx[ii, :], :]
        theta_mean = np.mean(tensor_normalize(theta_use, axis=1), axis=0)
        shape_sum  = np.sqrt(np.sum(theta_use[0, :]) * np.sum(theta_use[1, :]))
        theta_shapes2[ii, :] = theta_mean * shape_sum

    return np.append(theta_shapes, theta_shapes2, axis=0)


def add_doublet_GT(GT_prob):
    """
    Add doublet genotype by summarizing their probability:
    New GT has five categories: 0, 1, 2, 1.5, 2.5
    TODO: New GT has six categories: 0, 1, 2, 0_1, 0_2, 1_2
    """
    combn_iter = itertools.combinations(range(GT_prob.shape[1]), 2)
    gt_idx = np.array([x for x in combn_iter]) # GT combination
    g_idx1 = gt_idx[:, 0]
    g_idx2 = gt_idx[:, 1]

    combn_iter = itertools.combinations(range(GT_prob.shape[2]), 2)
    sp_idx = np.array([x for x in combn_iter]) # sample combination
    s_idx1 = sp_idx[:, 0]
    s_idx2 = sp_idx[:, 1]
    
    ## GT_prob has three genotypes: 0, 1, 2;
    n_gt = GT_prob.shape[1]
    GT_prob2 = np.zeros((GT_prob.shape[0], n_gt + gt_idx.shape[0], 
                         sp_idx.shape[0]))

    GT_prob2[:, :n_gt, :] = (GT_prob[:, :, s_idx1] * 
                             GT_prob[:, :, s_idx2])
    GT_prob2[:, n_gt:, :] = (GT_prob[:, :, s_idx1][:, g_idx1, :] * 
                             GT_prob[:, :, s_idx2][:, g_idx2, :] +
                             GT_prob[:, :, s_idx1][:, g_idx2, :] * 
                             GT_prob[:, :, s_idx2][:, g_idx1, :])
    
    GT_prob2 = tensor_normalize(GT_prob2, axis=1)
    GT_prob1 = np.append(GT_prob, np.zeros((GT_prob.shape[0], gt_idx.shape[0], 
                                            GT_prob.shape[2])), axis=1)
    return np.append(GT_prob1, GT_prob2, axis=2)
    

def greed_match(X, Z, axis=2):
    """
    Match Z to X by minimize the difference, 
    hence Z[:, :, axis] is best aligned to X
    """
    diff_mat = np.zeros((X.shape[axis], Z.shape[axis]))
    for i in range(X.shape[axis]):
        for j in range(Z.shape[axis]):
            diff_mat[i, j] = np.mean(np.abs(X[:, :, i] - Z[:, :, j]))
            
    diff_copy = diff_mat.copy()
    idx_out = -1 * np.ones(X.shape[axis], int)
    while (-1 in idx_out):
        idx_i = np.argmin(diff_copy) // diff_copy.shape[1]
        idx_j = np.argmin(diff_copy) % diff_copy.shape[1]
        idx_out[idx_i] = idx_j
        # print(idx_i, idx_j, idx_out)

        diff_copy[idx_i, :] = np.max(diff_mat) + 1
        diff_copy[:, idx_j] = np.max(diff_mat) + 1
        
    return idx_out


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
        _GT_diff = np.zeros((GT_prob.shape[2], GT_prob.shape[2]))
        for i in range(GT_prob.shape[2]):
            for j in range(GT_prob.shape[2]):
                _GT_diff[i, j] = np.mean(np.abs(GT_prob[:, :, i] - 
                                                GT_prob[:, :, j]))

        _donor_idx = [np.argmax(_donor_cnt)]
        _donor_left = np.delete(np.arange(GT_prob.shape[2]), _donor_idx)
        _GT_diff = np.delete(_GT_diff, _donor_idx, axis=1)
        while len(_donor_idx) < _GT_diff.shape[0]:
            # _idx = np.argmax(np.sum(_GT_diff[_donor_idx, :], axis=0))
            _idx = np.argmax(np.min(_GT_diff[_donor_idx, :], axis=0))
            _donor_idx.append(_donor_left[_idx])
            _donor_left = np.delete(_donor_left, _idx)
            _GT_diff = np.delete(_GT_diff, _idx, axis=1)

    print("\t".join(["donor%d" %x for x in _donor_idx]))
    print("\t".join(["%.0f" %_donor_cnt[x] for x in _donor_idx]))

    ID_prob_out = ID_prob[:, _donor_idx[:n_donor]]
    ID_prob_out[ID_prob_out < 10**-10] = 10**-10

    return ID_prob_out
