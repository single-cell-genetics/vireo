#Wrapper function for detecting useful mitochondrial variants

import os
from os import path
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
    def __init__(self, AD, DP, variant_names = None):
        #initiate object with AD/DP sparse matrices
        #check if AD and DP have same length first
        self.ad = AD.toarray()
        self.dp = DP.toarray()

        if len(self.ad) != len(self.dp):
            print('AD and DP length do not match!')
        else:
            print('%.2f variants detected' %(len(self.ad)))

        if variant_names is not None:
            #sanity check for length of variant names
            if len(variant_names) != len(self.ad):
                print('No. of variant names does not match length of AD!')
            else:
                self.variants = variant_names
        else:
            self.variants = None

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

        return len(_a), pval, params1, params2, model1.losses[-1], model2.losses[-1]

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

        return len(_a), p_val, params1, params2, model1.losses[-1], model2.losses[-1]
    
    def _deltaBIC(self, _a, _d, fix_seed=None, beta_mode=False):
        #input ad dp arrays, output params, BICs, delta BIC        
        if fix_seed is not None:
            np.random.seed(fix_seed)
        
        model1 = MixtureBinomial(n_components = 1, tor=1e-20)
        params1 = model1.fit((_a, _d), max_iters=500, early_stop=True)

        if beta_mode is False:
            model2 = MixtureBinomial(n_components = 2,tor=1e-20)
            params2 = model2.fit((_a, _d), max_iters=500, early_stop=True)
        else:
            model2 = MixtureBetaBinomial(n_components = 1, max_m_step_iter=3000,tor=1e-20, n_init_searches=100)
            params2 = model2.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)

        delta_BIC = model1.model_scores["BIC"] - model2.model_scores["BIC"]

        print("Cells qualified:%.2f\tmodel1 BIC:%.2f\tmodel2 BIC:%.2f\t deltaBIC:%.2f" %(len(_a),model2.model_scores["BIC"],model2.model_scores["BIC"],delta_BIC))

        return len(_a), delta_BIC, params1, params2, model1.model_scores["BIC"], model2.model_scores["BIC"]

    def fit_deltaBIC(self, nproc=30, minDP=10, minAD=1, beta_mode=False, export_csv = True):
        #here we fit and choose model based on deltaBIC
        print('CPUs used:', nproc)
        pool = mp.Pool(processes=nproc)
        results = []
        t0=time.time()

        print("Initializing fit(mode: deltaBIC) on %.2f variants..." %(len(self.ad)))

        for i in range(len(self.ad)):
            inputs = []
            idx = self.dp[i,:] >= minDP
            ad_idx = self.ad[i,:] >= minAD
            if any(idx) is True and any(ad_idx) is True:
                inputs.append([self.ad[i,idx], self.dp[i,idx], beta_mode])
                results.append(pool.starmap_async(self._deltaBIC, inputs))
            else:
                results.append(None)

        pool.close()
        pool.join()

        #num cells, deltaBIC, params1, params2, model1BIC, model2BIC
        self.output_list = [[] for i in range(6)]

        for res in results:
            if res is not None:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(res.get()[0][i])

            else:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(0)

        t1 = time.time()
        print("deltaBIC was calculated for %.2f variants and took:%.2f minutes" %(int(len(self.ad)),(t1-t0)/60))

        self.df = pd.DataFrame(data=self.output_list)
        self.df = self.df.transpose()
        self.df.columns = ['num_cells','deltaBIC', 'params1', 'params2', 'model1BIC', 'model2BIC']

        if export_csv is True:
            self.df.to_csv('BIC_params.csv', index=False)

        #return df of all metrics
        return self.df

    def fit_logLik(self, nproc=30, minDP=5, minAD=1, beta_mode=False, export_csv=True):
        #here we fit and choose model using Likelihood ratio test + neg log pval
        print('CPUs used:', nproc)
        pool = mp.Pool(processes=nproc)
        results = []
        t0 = time.time()
        print("Initializing fit(mode: LR) on %.2f variants..." %(len(self.ad)))

        if beta_mode is False:
            func = self._binomMixture
        else:
            func = self._betabinomMixture

        for i in range(len(self.ad)):
            inputs = []
            idx = self.dp[i,:] >= minDP
            ad_idx = self.ad[i,:] >= minAD
            if any(idx) is True and any(ad_idx) is True:
                inputs.append([self.ad[i,idx], self.dp[i,idx]])
                results.append(pool.starmap_async(func, inputs))
            else:
                results.append(None)

        pool.close()
        pool.join()

        #num cells, pval, params1, params2, model1 logLik, model2 logLik
        self.output_list = [[] for i in range(6)]

        for res in results:
            if res is not None:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(res.get()[0][i])

            else:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(0)

        t1 = time.time()
        print("deltaBIC was calculated for %.2f variants and took:%.2f minutes" %(int(len(self.ad)),(t1-t0)/60))

        self.df = pd.DataFrame(data=self.output_list)
        self.df = self.df.transpose()
        self.df.columns = ['num_cells','p_value', 'params1', 'params2', 'model1_logLik', 'model2_logLik']

        if export_csv is True:
            self.df.to_csv('pval_params.csv', index=False)

        #return df of all metrics
        return self.df

    def filter(self, by, threshold=500, export_heatmap=True, export_mtx=True, out_dir = None):
        #by should be a string

        if self.df is None:
            print('fitted model not found! Have you run fit_deltaBIC/fit_logLik yet?')
        else:
            if by == 'p_value':
                idx = np.where(np.array([-np.log(i) for i in self.df.p_value]) >= threshold)
            elif by == 'deltaBIC':
                idx = np.where(np.array(self.df.deltaBIC) >= threshold)

            best_ad = self.ad[idx]
            best_dp = self.dp[idx]
            print('Number of variants passing threshold:%.2f' %(len(best_ad)))

        fname = by + '_' + str(threshold) + '_'

        if self.variants is not None:
            best_vars = np.array(self.variants)[idx]

        if out_dir is not None:
            if path.exists(out_dir) is not True:
                try:
                    os.mkdir(out_dir)
                except:
                    print("Can't make directory, do you have permission?")
            else:
                print('Out directory already exists, overwriting content inside')
                
    
        if export_heatmap is True:
            af = best_ad/best_dp
            #af = af.fillna(0)
            fig, ax = plt.subplots(figsize=(15,10))
            plt.title("Allele frequency of top variants")
            plt.style.use('seaborn-dark')
            sns.heatmap(af, cmap='terrain_r')
            plt.savefig(out_dir + '/' + fname + 'top variants heatmap.pdf')

        #export ad dp mtx out for vireo
        if export_mtx is True:
            mmwrite(out_dir + '/' + fname + 'passed_ad.mtx', sparse.csr_matrix(best_ad))
            mmwrite(out_dir + '/' + fname + 'passed_dp.mtx', sparse.csr_matrix(best_dp))

        return best_ad, best_dp

if __name__ == '__main__':
    test_ad = mmread("C:/Users/aaronkwc/Downloads/dedup_kim_cellSNP.tag.AD.mtx")
    test_dp = mmread("C:/Users/aaronkwc/Downloads/dedup_kim_cellSNP.tag.DP.mtx")

    mdphd = MitoMut(AD = test_ad, DP = test_dp, variant_names = None)
    df = mdphd.fit_deltaBIC(nproc=30, beta_mode=False)
    print(df)
    final = mdphd.filter(by='deltaBIC', threshold=500, out_dir='kim_dataset_bic')
