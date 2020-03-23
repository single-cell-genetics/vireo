# Wrap function for Vireo model
# Author: Yuanhua Huang
# Date: 22/03/2020

import sys
import numpy as np
from .vireo_base import greed_match, donor_select
from .vireo_model import Vireo


def vireo_wrap(AD, DP, GT_prior=None, n_donor=None, learn_GT=True, n_init=20, 
    random_seed=None, check_doublet=True, max_iter_init=20, delay_fit_theta=3,
    n_extra_donor=0, extra_donor_mode="distance", **kwargs):
    """
    A wrap function to run vireo with multiple initializations
    """
    if learn_GT == False and n_extra_donor > 0:
        print("Searching from extra donors only works with learn_GT")
        n_extra_donor = 0
        
    # note learn_GT is false for mode 2 and 5 only (set before)
    if n_donor is None:
        if GT_prior is None:
            print("[vireo] Error: requiring n_donor or GT_prior.")
            sys.exit()
        else:
            n_donor = GT_prior.shape[1]

    if learn_GT is False and n_init > 1:
        print("GT is fixed, so use a single initialization")
        n_init = 1

    ## Setting random seed for initialization
    if random_seed is not None:
        np.random.seed(random_seed)

    GT_prior_use = None
    n_donor_use = int(n_donor + n_extra_donor)
    if GT_prior is not None and n_donor_use == GT_prior.shape[1]:
        GT_prior_use = GT_prior.copy()
    elif GT_prior is not None and n_donor_use < GT_prior.shape[1]:
        GT_prior_use = GT_prior.copy()
        n_donor_use = GT_prior.shape[1]

    _models_all = []
    for im in range(n_init):
        _modelCA = Vireo(n_var=AD.shape[0], n_cell=AD.shape[1], 
                         n_donor=n_donor_use, learn_GT=learn_GT, 
                         GT_prob_init=GT_prior_use, **kwargs)
        _modelCA.set_prior(GT_prior=GT_prior_use)
        _models_all.append(_modelCA)

    ## Fitting the models
    for im in range(n_init):
        # _models_all[im].fit(AD, DP, min_iter=20, verbose=False)
        _models_all[im].fit(AD, DP, min_iter=5, max_iter=max_iter_init, 
            delay_fit_theta=delay_fit_theta, verbose=False)

    ## select the model with best initialization
    elbo_all = np.array([x.ELBO_[-1] for x in _models_all])
    _idx = np.argmax(elbo_all)
    modelCA = _models_all[_idx]
    if n_extra_donor == 0:
        modelCA.fit(AD, DP, min_iter=5, verbose=False)
    else:
        _ID_prob = donor_select(modelCA.GT_prob, modelCA.ID_prob, n_donor, 
                                mode=extra_donor_mode)
        modelCA = Vireo(n_var=AD.shape[0], n_cell=AD.shape[1], 
                        n_donor=n_donor, learn_GT=learn_GT, 
                        GT_prob_init=GT_prior_use, ID_prob_init=_ID_prob,
                        beta_mu_init=modelCA.beta_mu, 
                        beta_sum_init=modelCA.beta_sum, **kwargs)
        modelCA.set_prior(GT_prior=GT_prior_use)
        modelCA.fit(AD, DP, min_iter=5, delay_fit_theta=delay_fit_theta, 
            verbose=False)

    print("[vireo] lower bound ranges [%.1f, %.1f, %.1f]" 
          %(np.min(elbo_all), np.median(elbo_all), np.max(elbo_all)))

    ## Run Vireo again with updateing genotype
    if GT_prior is not None and n_donor < GT_prior.shape[1]:
        _donor_cnt = np.sum(modelCA.ID_prob, axis=0)
        _donor_idx = np.argsort(_donor_cnt)[::-1]
        GT_prior_use = GT_prior[:, _donor_idx[:n_donor], :]

        modelCA = Vireo(n_var=AD.shape[0], n_cell=AD.shape[1], 
                        n_donor=n_donor, learn_GT=False, 
                        GT_prob_init=GT_prior_use, **kwargs)
        modelCA.fit(AD, DP, min_iter=20, verbose=False)

    elif GT_prior is not None and n_donor > GT_prior.shape[1]:
        GT_prior_use = modelCA.GT_prob.copy()
        idx = greed_match(GT_prior, GT_prior_use)
        GT_prior_use[:, idx, :] = GT_prior
        _idx_order = np.append(idx, np.delete(np.arange(n_donor), idx))
        GT_prior_use = GT_prior_use[:, _idx_order, :]
        ID_prob_use = modelCA.ID_prob[:, _idx_order]

        modelCA = Vireo(n_var=AD.shape[0], n_cell=AD.shape[1], 
                        n_donor=n_donor, learn_GT=learn_GT,
                        ID_prob_init=ID_prob_use,
                        beta_mu_init=modelCA.beta_mu, 
                        beta_sum_init=modelCA.beta_sum,
                        GT_prob_init=GT_prior_use, **kwargs)
        modelCA.set_prior(GT_prior = GT_prior_use)
        modelCA.fit(AD, DP, min_iter=20, verbose=False)
    
    ## print the beta parameters
    print("[vireo] allelic rate mean and concentrations:")
    print(np.round(modelCA.beta_mu, 3))
    print(np.round(modelCA.beta_sum, 1))

    ## Summarise donor size
    print("[vireo] donor size before removing doublets:")
    _donor_cnt = np.sum(modelCA.ID_prob, axis=0)
    print("\t".join(["donor%d" %x for x in range(len(_donor_cnt))]))
    print("\t".join(["%.0f" %x for x in _donor_cnt]))
    
    
    ## Predict doublets
    if check_doublet:
        doublet_prob, ID_prob = modelCA.predict_doublet(AD, DP)
    else:
        ID_prob = modelCA.ID_prob
        doublet_prob = np.zeros((AD.shape[1], AD.shape[1] * (AD.shape[1] - 1) / 2))

    theta_shapes = np.append(modelCA.beta_mu * modelCA.beta_sum, 
                             (1 - modelCA.beta_mu) * modelCA.beta_sum, axis=0)
    RV = {}
    RV['ID_prob'] = ID_prob
    RV['GT_prob'] = modelCA.GT_prob
    RV['doublet_prob'] = doublet_prob
    RV['theta_shapes'] = theta_shapes
    RV['theta_mean'] = modelCA.beta_mu
    RV['theta_sum'] = modelCA.beta_sum
    RV['LB_list'] = elbo_all
    RV['LB_doublet'] = modelCA.ELBO_[-1]
    return RV
