#Wrapper function for detecting useful mitochondrial variants

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread
from scipy.io import mmwrite
from scipy import sparse
from scipy.stats import betabinom, bernoulli, binom

import bbmix
from bbmix.models import MixtureBinomial
from bbmix.models import MixtureBetaBinomial

import multiprocessing as mp

class MitoMut():
    def __init__(self, AD, DP, variant_names):
        #initiate object with AD/DP sparse matrices
        self.ad = AD.toarray()
        self.dp = DP.toarray()
        self.variants = variant_names

    def _betabinomMixture(self, _a, _d, fix_seed=False):
        #input ad dp arrays, output pval
        model1 = MixtureBetaBinomial(n_components = 1, max_m_step_iter=3000,tor=1e-20, n_init_searches=100)
        model2 = MixtureBetaBinomial(n_components = 2, max_m_step_iter=3000,tor=1e-20, n_init_searches=500)

        if fix_seed is True:
            np.random.seed(42)

        params1 = model1.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)
        params2 = model2.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)
        p_val = bbmix.models.LR_test(model1.losses[-1] - model2.losses[-1], df = 3)
        print("Cells qualified:%.2f\tmodel1:%.2f\tmodel2:%.2f\tp value:%.2f" %(len(_a),model1.losses[-1],model2.losses[-1],p_val))

        return p_val

    def _binomMixture(self, _a, _d, fix_seed=False):
        #input ad dp arrays, output pval
        model1 = MixtureBinomial(n_components = 1, tor=1e-20)
        model2 = MixtureBinomial(n_components = 2,tor=1e-20)

        if fix_seed is True:
            np.random.seed(42)

        params1 = model1.fit((_a, _d), max_iters=500, early_stop=True)
        params2 = model2.fit((_a, _d), max_iters=500, early_stop=True)
        p_val = bbmix.models.LR_test(model1.losses[-1] - model2.losses[-1], df = 2)
        print("Cells qualified:%.2f\tmodel1:%.2f\tmodel2:%.2f\tp value:%.2f" %(len(_a),model1.losses[-1],model2.losses[-1],p_val))

        return p_val

    def fitBinom(self, nproc=30, minDP=5, minAD=1):

        pool = mp.Pool(processes=nproc)
        results = []
        t0 = time.time()
        print("Initializing binomial mixture fit on %.2f variants..." %(len(self.ad)))

        for i in range(len(self.ad)):
            inputs = []
            idx = self.dp[i,:] >= minDP
            ad_idx = self.ad[i,:] >= minAD
            if any(idx) is True and any(ad_idx) is True:
                inputs.append([self.ad[i,idx], self.dp[i,idx]])
                results.append(pool.starmap_async(self._binomMixture, inputs))
            else:
                results.append(None)

        pool.close()
        pool.join()

        self.pvalue_b = []
        for res in results:
            if res is not None:
                self.pvalue_b.append(res.get())
            else:
                self.pvalue_b.append(1)

        t1 = time.time()
        print("Binomial mixture model processed %.2f variants and took:%.2f minutes" %(int(len(self.ad)),(t1-t0)/60))

        #return a list of pvalues
        return self.pvalue_b

    def filterBinom(self, threshold=100):
        log_pv = [-np.log(res) for res in self.pvalue_b]
        idx = np.where(np.array(log_pv) >= threshold)[0]
        
        self.top_ad = self.ad[idx]
        self.top_dp = self.dp[idx]

        if self.variants is not None:
            self.topvars = np.array(self.variants)[idx]
        else:
            self.topvars = None

        print(str(log_pv))

        print("%.2f variants passed the threshold!" %(len(self.top_ad)))

        return self.top_ad, self.top_dp, self.topvars 

    def fitBetaBinom(self, nproc=30, minDP=5, minAD=1):

        pool = mp.Pool(processes=nproc)
        results = []
        t0 = time.time()
        print("Initializing beta binomial mixture fit on %.2f variants..." %(len(self.top_ad)))

        for i in range(len(self.top_ad)):
            inputs = []
            idx = self.top_dp[i,:] >= minDP
            inputs.append([self.top_ad[i,idx], self.top_dp[i,idx]])
            results.append(pool.starmap_async(self._betabinomMixture, inputs))

        pool.close()
        pool.join()

        self.pvalue_bb = []
        for res in results:
            if res is not None:
                self.pvalue_bb.append(res.get())
            else:
                self.pvalue_bb.append(1)

        t1 = time.time()
        print("BetaBinomial model looped through %.2f variants and took:%.2f minutes" %(int(len(self.top_ad)),(t1-t0)/60))

        #return a list of pvalues
        return self.pvalue_bb

    def filterBetaBinom(self, threshold=5):
        log_pv = [-np.log(res) for res in self.pvalue_bb]
        idx = np.where(np.array(log_pv) >= threshold)[0]
        
        self.best_ad = self.top_ad[idx]
        self.best_dp = self.top_dp[idx]

        if self.topvars is not None:
            self.bestvars = np.array(self.topvars)[idx]
        else:
            self.bestvars = None

        return self.best_ad, self.best_dp, self.bestvars

"""
if __name__ == '__main__':
    test_ad = mmread("data/mitoDNA/cellSNP.tag.AD.mtx")
    test_dp = mmread("data/mitoDNA/cellSNP.tag.DP.mtx")

    mdphd = MitoMut(AD = test_ad, DP = test_dp, variant_names = None)
    mdphd.fitBinom(nproc = 15)
    res_bin = mdphd.filterBinom()
    mdphd.fitBetaBinom(nproc = 30)
    res_bb = mdphd.filterBetaBinom()
"""