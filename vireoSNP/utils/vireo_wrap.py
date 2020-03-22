import numpy as np
from .vireo_base import normalize, greed_match, donor_select
from .vireo_model import vireo_core


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
    if GT_prior is None or n_donor > GT_prior.shape[1]:
        n_donor_run1 = int(n_donor + n_extra_donor)
        print("[vireo] warm-up: %d random initializations for %d clusters..." 
              %(n_init, n_donor_run1))

        ID_prob_list = []
        for i in range(n_init):
            _ID_prob = np.random.rand(AD.shape[1], n_donor_run1)
            ID_prob_list.append(normalize(_ID_prob))

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
    if GT_prior is not None and n_donor < GT_prior.shape[1]:
        print("[vireo] pre-RUN: finding %d from %d donors with GT ..." 
              %(n_donor, GT_prior.shape[2]))
        res1 = vireo_core(AD, DP, GT_prior = GT_prior, n_donor=None, 
                          ID_prob_init=_ID_prob, check_doublet=False, **kwargs)
        _donor_cnt = np.sum(res1['ID_prob'], axis=0)
        _donor_idx = np.argsort(_donor_cnt)[::-1]
        GT_prior_use = GT_prior[:, _donor_idx[:n_donor], :]

        print("\t".join(["donor%d" %x for x in _donor_idx]))
        print("\t".join(["%.0f" %_donor_cnt[x] for x in _donor_idx]))

    elif GT_prior is not None and n_donor > GT_prior.shape[1]:
        print("[vireo] pre-RUN: finding %d from %d donors without GT ..." 
              %(n_donor - GT_prior.shape[1], n_donor))
        res1 = vireo_core(AD, DP, GT_prior=None, n_donor=n_donor, 
                          ID_prob_init=_ID_prob, check_doublet=check_doublet, 
                          **kwargs)
        GT_prior_use = res1['GT_prob']
        idx = greed_match(GT_prior, GT_prior_use)
        GT_prior_use[:, idx, :] = GT_prior
        _idx_order = np.append(idx, np.delete(np.arange(n_donor), idx))
        GT_prior_use = GT_prior_use[:, _idx_order, :]

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
    np.set_printoptions(formatter={'float': lambda x: format(x, '.2f')})
    print(res1['theta_shapes'])

    return res1
