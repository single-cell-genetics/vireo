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
    def __init__(self, AD, DP, variant_names = None, dataset_name = None):
        #initiate object with AD/DP sparse matrices
        #check if AD and DP have same length first
        self.ad = AD.toarray()
        self.dp = DP.toarray()

        if len(self.ad) != len(self.dp):
            print('AD and DP length do not match!')
        else:
            print(str(len(self.ad)) + ' variants detected')

        if variant_names is not None:
            #sanity check for length of variant names
            if len(variant_names) != len(self.ad):
                print('No. of variant names does not match length of AD!')
            else:
                self.variants = variant_names
                print("variant names detected") 
            
        else:
            self.variants = None

        if dataset_name is not None:
            self.dataset = dataset_name

    def _betabinomMixture(self, _a, _d, fix_seed=False):
        #basic staistics
        #Total DP across all cells
        total_DP = np.sum(_d)
        #Median DP across all cells
        median_DP = np.median(_d)
        #Total AD across all cells
        total_AD = np.sum(_a)
        #Median AD across all cells
        median_AD = np.median(_a)
        #How many cells have this variant?
        non_zero = np.count_nonzero(_a)

        #input ad dp arrays, output pval
        model1 = MixtureBetaBinomial(n_components = 1, max_m_step_iter=3000,tor=1e-20, n_init_searches=100)
        model2 = MixtureBetaBinomial(n_components = 2, max_m_step_iter=3000,tor=1e-20, n_init_searches=500)

        if fix_seed is True:
            np.random.seed(42)

        params1 = model1.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)
        params2 = model2.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)
        p_val = bbmix.models.LR_test(model1.losses[-1] - model2.losses[-1], df = 3)
        print("Cells qualified: " + str(len(_a)) + "\tmodel1:%.2f\tmodel2:%.2f\tp value:%.2f" %(model1.losses[-1],model2.losses[-1],p_val))

        return len(_a), p_val, params1, params2, model1.losses[-1], model2.losses[-1], non_zero, total_DP, median_DP, total_AD, median_AD

    def _binomMixture(self, _a, _d, fix_seed=False):
        #basic staistics
        #Total DP across all cells
        total_DP = np.sum(_d)
        #Median DP across all cells
        median_DP = np.median(_d)
        #Total AD across all cells
        total_AD = np.sum(_a)
        #Median AD across all cells
        median_AD = np.median(_a)
        #How many cells have this variant?
        non_zero = np.count_nonzero(_a)

        #input ad dp arrays, output pval
        model1 = MixtureBinomial(n_components = 1, tor=1e-20)
        model2 = MixtureBinomial(n_components = 2,tor=1e-20)

        if fix_seed is True:
            np.random.seed(42)

        params1 = model1.fit((_a, _d), max_iters=500, early_stop=True)
        params2 = model2.fit((_a, _d), max_iters=500, early_stop=True)
        p_val = bbmix.models.LR_test(model1.losses[-1] - model2.losses[-1], df = 2)
        print("Cells qualified: " + str(len(_a)) + "\tmodel1:%.2f\tmodel2:%.2f\tp value:%.2f" %(model1.losses[-1],model2.losses[-1],p_val))

        return len(_a), p_val, params1, params2, model1.losses[-1], model2.losses[-1], non_zero, total_DP, median_DP, total_AD, median_AD
    
    def _deltaBIC(self, _a, _d, fix_seed=None, beta_mode=False, sim=False):
        #input ad dp arrays, output params, BICs, delta BIC        
        if fix_seed is not None:
            np.random.seed(fix_seed)
        
        #basic staistics
        #Total DP across all cells
        total_DP = np.sum(_d)
        #Median DP across all cells
        median_DP = np.median(_d)
        #Total AD across all cells
        total_AD = np.sum(_a)
        #Median AD across all cells
        median_AD = np.median(_a)
        #How many cells have this variant?
        non_zero = np.count_nonzero(_a)

        model1 = MixtureBinomial(n_components = 1, tor=1e-20)
        params1 = model1.fit((_a, _d), max_iters=500, early_stop=True)

        if beta_mode is False:
            model2 = MixtureBinomial(n_components = 2,tor=1e-20)
            if sim is False:
                params2 = model2.fit((_a, _d), max_iters=500, early_stop=True)
            else:
                a_sim = model1.sample(_d)
                params1 = model1.fit((a_sim, _d), max_iters=500, early_stop=True)
                params2 = model2.fit((a_sim, _d), max_iters=500, early_stop=True)
        else:
            model2 = MixtureBetaBinomial(n_components = 1, max_m_step_iter=3000,tor=1e-20, n_init_searches=100)
            params2 = model2.fit((_a, _d), max_iters=3000, init_method="mixbin", early_stop=False, n_tolerance=10)

        delta_BIC = model1.model_scores["BIC"] - model2.model_scores["BIC"]

        p = params2[0] , params2[1]
        pi = params2[2], params2[3]
        fraction_b_allele = np.min(np.array(p)) * np.array(pi)[np.argmin(np.array(p))]

        if np.max(np.array(pi)) < 0.95 and np.min(np.array(p)) < 0.05 and np.max(np.array(p)) > 0.1:
            new_mutation = True
            as_mutation = False
        elif np.min(np.array(p)) > 0.1 and np.min(np.array(pi)) > 0.15:
            as_mutation = True
            new_mutation = False
        else:
            new_mutation, as_mutation = False, False

        print("Cells qualified: " + str(len(_a)) + "\tmodel1 BIC:%.2f\tmodel2 BIC:%.2f\t deltaBIC:%.2f" %(model2.model_scores["BIC"],model2.model_scores["BIC"],delta_BIC))

        return len(_a), delta_BIC, params1, params2, model1.model_scores["BIC"], model2.model_scores["BIC"], non_zero, total_DP, median_DP, total_AD, median_AD, new_mutation, as_mutation, fraction_b_allele

    def _check_outdir_exist(self, out_dir):
        if path.exists(out_dir) is not True:
            try:
                os.mkdir(out_dir)
                return True
            except:
                print("Can't make directory, do you have permission?")
                return False
        else:
            print('Out directory already exists, overwriting content inside')
            return True

    def fit_deltaBIC(self, out_dir, nproc=30, minDP=10, minAD=1, beta_mode=False, export_csv = True, sim=False):
        #here we fit and choose model based on deltaBIC
        print('CPUs used:', nproc)
        pool = mp.Pool(processes=nproc)
        results = []
        t0=time.time()

        print("Initializing fit(mode: deltaBIC) on " + str(len(self.ad)) + " variants...")

        for i in range(len(self.ad)):
            inputs = []
            idx = self.dp[i,:] >= minDP
            ad_idx = self.ad[i,:] >= minAD
            if any(idx) is True and any(ad_idx) is True:
                inputs.append([self.ad[i,idx], self.dp[i,idx], beta_mode, sim])
                results.append(pool.starmap_async(self._deltaBIC, inputs))
            else:
                results.append(None)

        pool.close()
        pool.join()

        #num cells, deltaBIC, params1, params2, model1BIC, model2BIC
        self.output_list = [[] for i in range(14)]

        for res in results:
            if res is not None:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(res.get()[0][i])

            else:
                for i in range(len(self.output_list)):
                    self.output_list[i].append(0)

        t1 = time.time()
        print("deltaBIC was calculated for " + str(len(self.ad)) + " variants and took:%.2f minutes" %((t1-t0)/60))

        self.df = pd.DataFrame(data=self.output_list)
        self.df = self.df.transpose()
        self.df.columns = ['num_cells','deltaBIC', 'params1', 'params2', 'model1BIC', 'model2BIC', 'num_cells_nonzero_AD', 'total_DP', 'median_DP', 'total_AD', 'median_AD', 'new_mutation', 'as_mutation', 'fraction_b_allele']

        if self.variants is not None:
            self.df = pd.concat([pd.Series(self.variants), self.df], axis=1)

        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        if export_csv is True:
            if self._check_outdir_exist(out_dir) is True:
                self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)
            else:
                self.sorted_df.to_csv('BIC_params.csv', index=False)

        self.df.to_csv(out_dir + '/debug_unsorted_BIC_params.csv', index=False)
        #return df of all metrics
        return self.df

    def fit_logLik(self, out_dir, nproc=30, minDP=5, minAD=1, beta_mode=False, export_csv=True):
        #here we fit and choose model using Likelihood ratio test + neg log pval
        print('CPUs used:', nproc)
        pool = mp.Pool(processes=nproc)
        results = []
        t0 = time.time()
        print("Initializing fit(mode: LR) on" + len(self.ad) + "variants...")

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
        print("Log likelihood was calculated for" + str(len(self.ad)) + " variants and took:%.2f minutes" %((t1-t0)/60))

        self.df = pd.DataFrame(data=self.output_list)
        self.df = self.df.transpose()
        self.df.columns = ['num_cells','p_value', 'params1', 'params2', 'model1_logLik', 'model2_logLik']
        
        if self.variants is not None:
            self.df = pd.concat([pd.Series(self.variants), self.df], axis=1)

        self.sorted_df = self.df.sort_values(by=['p_value'])

        if export_csv is True:
            if self._check_outdir_exist(out_dir) is True:
                self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)
            else:
                self.sorted_df.to_csv('pval_params.csv', index=False)

        #return df of all metrics
        return self.df

    def filter(self, by, threshold=500, export_heatmap=True, export_mtx=True, out_dir=None):
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
            print('Number of variants passing threshold: '  + str(len(best_ad)))

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
            if self.variants is not None:
                sns.heatmap(af, cmap='terrain_r', yticklabels=best_vars)
            else:
                sns.heatmap(af, cmap='terrain_r')
            plt.savefig(out_dir + '/' + fname + 'top variants heatmap.pdf')

        #export ad dp mtx out for vireo
        if export_mtx is True:
            mmwrite(out_dir + '/' + fname + 'passed_ad.mtx', sparse.csr_matrix(best_ad))
            mmwrite(out_dir + '/' + fname + 'passed_dp.mtx', sparse.csr_matrix(best_dp))

        return best_ad, best_dp

    def rankBIC(self, top, export_heatmap=True, export_mtx=True, out_dir=None):
        #test function for ranking variants instead of using a cutoff
        #top_rank = zip(*heapq.nlargest(top, enumerate(self.df.deltaBIC), key=operator.itemgetter(1)))[0]
        self.new_mut_df = self.df[(self.df.new_mutation == 'True')|(self.df.as_mutation == 'True')]
        #self.new_mut_df = self.df[self.df.new_mutation == 'True']
        #self.as_mut_df = self.df[self.df.as_mutation == 'True'] 
        #print(self.new_mut_df)
        #top_rank = np.argpartition(np.array(self.new_mut_df.deltaBIC), -top)[-top:]
        #new_top_rank = self.new_mut_df.sort_values(by='deltaBIC', ascending=False)[0:top].index
        #as_top_rank = self.as_mut_df.sort_values(by='deltaBIC', ascending=False)[0:top].index

        idx = self.new_mut_df[self.new_mut_df.deltaBIC >= 500].index
        print(len(idx))
        #print(idx)
        #idx2 = self.as_mut_df[self.as_mut_df.deltaBIC >= 500].index
        #print(idx2)
        #both_top_rank = idx.append(idx2)
        #print(both_top_rank)
        best_ad, best_dp = self.ad[idx], self.dp[idx]



        fname = 'top_' + str(top) + '_'

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
            if self.variants is not None:
                sns.heatmap(af, cmap='terrain_r', yticklabels=best_vars)
            else:
                sns.heatmap(af, cmap='terrain_r')
            plt.savefig(out_dir + '/' + fname + 'variants heatmap.pdf')

        #export ad dp mtx out for vireo
        if export_mtx is True:
            mmwrite(out_dir + '/' + fname + 'passed_ad.mtx', sparse.csr_matrix(best_ad))
            mmwrite(out_dir + '/' + fname + 'passed_dp.mtx', sparse.csr_matrix(best_dp))

        return best_ad, best_dp
    
    def read_df(self, file):
        self.df = pd.read_csv(file)
        return self.df

if __name__ == '__main__':
    import vireoSNP
    from vireoSNP import __version__
    
    #from vireoSNP.utils.io_utils import read_sparse_GeneINFO
    #from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb

    #test_ad = mmread("data/mitoDNA/cellSNP.tag.AD.mtx")
    #test_dp = mmread("data/mitoDNA/cellSNP.tag.DP.mtx")

    #cell_vcf = load_VCF("data/mitoDNA/kim_cellSNP.cells.vcf.gz", biallelic_only=True)
    #cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
    #for _key in ['samples', 'variants', 'FixedINFO', 'contigs', 'comments']:
        #cell_dat[_key] = cell_vcf[_key]
    
    cell_vcf = vireoSNP.load_VCF("data/kim/cellSNP.cells.vcf.gz", biallelic_only=True)
    cell_dat = vireoSNP.vcf.read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
    mdphd = MitoMut(AD = cell_dat['AD'], DP = cell_dat['DP'], 
                    variant_names = cell_vcf['variants'])

    #mdphd = MitoMut(AD = test_ad, DP = test_dp)
    df = mdphd.fit_deltaBIC(out_dir='data/kim/sim_out', nproc=15, sim=True)
    #final = mdphd.filter(by='deltaBIC', threshold = 500, out_dir = 'data/mitoDNA/mitoMutOUT')
    #df = mdphd.read_df('data/MGH26/out/debug_unsorted_BIC_params.csv')
    #final = mdphd.rankBIC(top=500, out_dir='data/MGH26/out')
    #final2 = mdphd.rankBIC(top=100, out_dir='data/mitoDNA/mitoMutOUT')
    #final3 = mdphd.rankBIC(top=50, out_dir='data/mitoDNA/mitoMutOUT')