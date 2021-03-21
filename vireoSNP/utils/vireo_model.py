import itertools
import numpy as np
from scipy.stats import entropy
from scipy.sparse import csc_matrix
from scipy.special import logsumexp, digamma, betaln
from .vireo_base import normalize, loglik_amplify, beta_entropy
from .vireo_base import get_binom_coeff, logbincoeff

__docformat__ = "restructuredtext en"

class Vireo():
    """Viroe model: Variational Inference for reconstruction of ensemble origin

    The prior can be set via set_prior() before fitting the model.

    Key properties
    --------------
    beta_mu: numpy array (1, n_GT) or (n_var, n_GT)
        Beta mean parameter of theta's posterior
    beta_sum: numpy array (1, n_GT) or (n_var, n_GT), same as beta_mu
        Beta concetration parameter of theta's posterior
    ID_prob: numpy array (n_cell, n_donor)
        Posterior cell assignment probability to each donor
    GT_prob: numpy array (n_var, n_donor, n_GT)
        Posterior genotype probability per variant per donor
    """
    def __init__(self, n_cell, n_var, n_donor, n_GT=3, learn_GT=True,  
        learn_theta=True, ASE_mode=False, fix_beta_sum=False, 
        beta_mu_init=None, beta_sum_init=None, ID_prob_init=None, 
        GT_prob_init=None):
        """Initialise Vireo model

        Note, multiple initializations are highly recomended to avoid local 
        optima.
        
        Parameters
        ----------
        n_cell : int. 
            Number of cells
        n_var : int. 
            Number of variants
        n_donor : int. 
            Number of donors
        n_GT : int. 
            Number of genotype categories
        learn_GT: bool. 
            Whether updating `GT_prob`; otherwise using the initial
        ASE_mode: bool. 
            Whether setting allelic ratio `theta` to be variant specific
        fix_beta_sum: bool. 
            Whether fixing the concetration parameter of theta's posterior
        beta_mu_init: numpy array (1, n_GT) or (n_var, n_GT)
            Initial value of beta_mu, the mean parameter of theta
        beta_sum_init: numpy array (1, n_GT) or (n_var, n_GT), same as beta_mu
            Initial value of beta_sum, the concetration parameter of theta
        ID_prob_init: numpy array (n_cell, n_donor)
            Initial value of ID_prob, cell assignment probability to each donor
        GT_prob_init: numpy array (n_var, n_donor, n_GT)
            Initial value of GT_prob, genotype probability per variant and donor
        """
        self.n_GT = n_GT
        self.n_var = n_var
        self.n_cell = n_cell
        self.n_donor = n_donor
        self.learn_GT = learn_GT
        self.ASE_mode = ASE_mode
        self.learn_theta = learn_theta
        self.fix_beta_sum = fix_beta_sum
        
        self.ELBO_ = np.zeros((0))
        
        # initial key parameters
        self.set_initial(beta_mu_init, beta_sum_init, ID_prob_init, GT_prob_init)
        
        # set hyper parameters for prior
        self.set_prior()

    def set_initial(self, beta_mu_init=None, beta_sum_init=None, 
        ID_prob_init=None, GT_prob_init=None):
        """Set initial values
        """
        theta_len = self.n_var if self.ASE_mode else 1

        if beta_mu_init is not None:
            self.beta_mu = beta_mu_init
        else:
            self.beta_mu = (np.ones((theta_len, self.n_GT)) * 
                np.linspace(0.01, 0.99, self.n_GT).reshape(1, -1))

        if beta_sum_init is not None:
            self.beta_sum = beta_sum_init
        else:
            self.beta_sum = np.ones((theta_len, self.n_GT)) * 50

        if ID_prob_init is not None:
            self.ID_prob = normalize(ID_prob_init, axis=1)
        else:
            self.ID_prob = normalize(np.random.rand(self.n_cell, self.n_donor))

        if GT_prob_init is not None:
            self.GT_prob = normalize(GT_prob_init)
        else:        
            _GT_val = np.random.rand(self.n_var, self.n_donor, self.n_GT)
            self.GT_prob = normalize(_GT_val)

    
    def set_prior(self, GT_prior=None, ID_prior=None, beta_mu_prior=None, 
        beta_sum_prior=None, min_GP=0.00001):
        """Set prior for key variables: theta, GT_prob and ID_prob.
        The priors are in the same shape as its according variables.

        min_GP: float. Minimun genotype probability in GT_prior.
        """
        if beta_mu_prior is None:
            beta_mu_prior = np.expand_dims(
                np.linspace(0.01, 0.99, self.beta_mu.shape[1]), axis=0)
        if beta_sum_prior is None:
            beta_sum_prior = np.ones(beta_mu_prior.shape) * 50.0
        self.theta_s1_prior = beta_mu_prior * beta_sum_prior
        self.theta_s2_prior = (1 - beta_mu_prior) * beta_sum_prior

        if ID_prior is not None:
            if len(ID_prior.shape) == 1:
                ID_prior = np.expand_dims(ID_prior, axis=0)
            self.ID_prior = ID_prior
        else:
            self.ID_prior = normalize(np.ones(self.ID_prob.shape))

        if GT_prior is not None:
            if len(GT_prior.shape) == 2:
                GT_prior = np.expand_dims(GT_prior, axis=0)
            GT_prior[GT_prior < min_GP] = min_GP
            GT_prior[GT_prior > 1 - min_GP] = 1 - min_GP
            GT_prior = normalize(GT_prior)
            self.GT_prior = GT_prior
        else:        
            self.GT_prior = normalize(np.ones(self.GT_prob.shape))

    @property
    def theta_s1(self):
        """Beta concetration1 parameter for theta posterior"""
        return self.beta_mu * self.beta_sum

    @property
    def theta_s2(self):
        """Beta concetration2 parameter for theta posterior"""
        return (1 - self.beta_mu) * self.beta_sum

    @property 
    def digamma1_(self):
        """Digamma of Beta concetration1 parameter"""
        return np.expand_dims(digamma(self.theta_s1), 1)

    @property 
    def digamma2_(self):
        """Digamma of Beta concetration2 parameter"""
        return np.expand_dims(digamma(self.theta_s2), 1)

    @property 
    def digammas_(self):
        """Digamma of Beta concetration summary parameter"""
        return np.expand_dims(digamma(self.theta_s1 + self.theta_s2), 1)


    def update_theta_size(self, AD, DP):
        """Coordinate ascent for updating theta posterior parameters
        """
        BD = DP - AD
        S1_gt = AD @ self.ID_prob  #(n_var, n_donor)
        S2_gt = BD @ self.ID_prob  #(n_var, n_donor)
        
        _theta_s1 = np.zeros(self.beta_mu.shape)
        _theta_s2 = np.zeros(self.beta_mu.shape)
        _theta_s1 += self.theta_s1_prior.copy()
        _theta_s2 += self.theta_s2_prior.copy()
        for ig in range(self.n_GT):
            _axis = 1 if self.ASE_mode else None
            _theta_s1[:, ig:(ig+1)] += np.sum(
                S1_gt * self.GT_prob[:, :, ig], axis=_axis, keepdims=True)
            _theta_s2[:, ig:(ig+1)] += np.sum(
                S2_gt * self.GT_prob[:, :, ig], axis=_axis, keepdims=True)
        
        self.beta_mu = _theta_s1 / (_theta_s1 + _theta_s2)
        if self.fix_beta_sum == False:
            self.beta_sum = _theta_s1 + _theta_s2

    def update_ID_prob(self, AD, DP):
        """Coordinate ascent for updating assignment probability
        """
        BD = DP - AD
        logLik_ID = np.zeros((AD.shape[1], self.n_donor))
        for ig in range(self.n_GT):
            S1 = AD.T @ (self.GT_prob[:, :, ig] * self.digamma1_[:, :, ig])
            S2 = BD.T @ (self.GT_prob[:, :, ig] * self.digamma2_[:, :, ig])
            SS = DP.T @ (self.GT_prob[:, :, ig] * self.digammas_[:, :, ig])
            logLik_ID += (S1 + S2 - SS)
        
        self.ID_prob = normalize(np.exp(loglik_amplify(
            logLik_ID + np.log(self.ID_prior))))
        
        return logLik_ID
                

    def update_GT_prob(self, AD, DP):
        """Coordinate ascent for updating genotype probability
        """
        S1_gt = AD @ self.ID_prob
        SS_gt = DP @ self.ID_prob
        S2_gt = SS_gt - S1_gt
        
        logLik_GT = np.zeros(self.GT_prior.shape)
        for ig in range(self.n_GT):
            logLik_GT[:, :, ig] = (
                S1_gt * self.digamma1_[:, :, ig] + 
                S2_gt * self.digamma2_[:, :, ig] - 
                SS_gt * self.digammas_[:, :, ig])

        self.GT_prob = normalize(np.exp(loglik_amplify(
            logLik_GT + np.log(self.GT_prior))))
        
        
    def get_ELBO(self, logLik_ID, AD=None, DP=None):
        """Calculating variational evidence lower bound with current parameters

        logLik_ID: numpy array (n_cell, n_donor), the output from update_ID_prob
        """
        if logLik_ID is None:
            BD = DP - AD
            logLik_ID = np.zeros((AD.shape[1], self.n_donor))
            for ig in range(self.n_GT):
                S1 = AD.T @ (self.GT_prob[:, :, ig] * self.digamma1_[:, :, ig])
                S2 = BD.T @ (self.GT_prob[:, :, ig] * self.digamma2_[:, :, ig])
                SS = DP.T @ (self.GT_prob[:, :, ig] * self.digammas_[:, :, ig])
                logLik_ID += (S1 + S2 - SS)

        LB_p = np.sum(logLik_ID * self.ID_prob)
        KL_ID = np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1))
        KL_GT = np.sum(entropy(self.GT_prob, self.GT_prior, axis=-1))
        KL_theta = beta_entropy(
            np.append(
                np.expand_dims(self.theta_s1, 1), 
                np.expand_dims(self.theta_s2, 1), axis = 1),
            np.append(
                np.expand_dims(self.theta_s1_prior, 1), 
                np.expand_dims(self.theta_s2_prior, 1), axis = 1))
        
        # print(LB_p, KL_ID, KL_GT, KL_theta)
        return LB_p - KL_ID - KL_GT - KL_theta


    def _fit_VB(self, AD, DP, max_iter=200, min_iter=5, epsilon_conv=1e-2,
        delay_fit_theta=0, verbose=True):
        """Fit Vireo model with coordinate ascent
        """
        ELBO = np.zeros(max_iter)
        numerical_minimal = 1e-6
        for it in range(max_iter):
            if self.learn_theta and it >= delay_fit_theta:
                self.update_theta_size(AD, DP)
            if self.learn_GT:
                self.update_GT_prob(AD, DP)

            _logLik_ID = self.update_ID_prob(AD, DP)
            ELBO[it] = self.get_ELBO(_logLik_ID) #+ _binom_coeff_log
            
            if it > min_iter:
                if ELBO[it] < ELBO[it - 1] - numerical_minimal:
                    if verbose:
                        print("Warning: Lower bound decreases!\n")
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif ELBO[it] - ELBO[it - 1] < epsilon_conv:
                    break
        
        return ELBO[:it]

    def fit(self, AD, DP, max_iter=200, min_iter=5, epsilon_conv=1e-2,
        delay_fit_theta=0, verbose=True, n_inits=50, nproc=1):
        """Fit Vireo model with coordinate ascent

        Parameters
        ----------
        AD : scipy.sparse.csc_matrix (n_var, n_cell)
            Sparse count matrix for alternative allele
        DP : scipy.sparse.csc_matrix (n_var, n_cell)
            Sparse count matrix for depths, alternative + refeerence alleles
        max_iter : int
            Maximum number of iterations
        min_iter :
            Minimum number of iterations
        epsilon_conv : float
            Threshold for detecting convergence
        delay_fit_theta : int
            Number of steps to delay updating theta. This can be very useful 
            for common genetics when there is good prior on allelic ratio.
        verbose : bool
            Whether print out log info
        """
        if type(DP) is np.ndarray and np.mean(DP > 0) < 0.3:
            print("Warning: input matrices is %.1f%% sparse, " 
                  %(100 - np.mean(DP > 0) * 100) +
                  "change to scipy.sparse.csc_matrix" )
            AD = csc_matrix(AD)
            DP = csc_matrix(DP)

        ELBO = self._fit_VB(AD, DP, max_iter, min_iter, epsilon_conv, 
                            delay_fit_theta, verbose)

        # _binom_coeff_log = np.sum(logbincoeff(DP, AD, is_sparse=True))
        # _binom_coeff_log = np.sum(get_binom_coeff(AD, DP))
        
        ELBO += np.sum(get_binom_coeff(AD, DP))

        self.ELBO_ = np.append(self.ELBO_, ELBO)
