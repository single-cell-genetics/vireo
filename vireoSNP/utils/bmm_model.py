import itertools
import numpy as np
from scipy.stats import entropy
from scipy.sparse import csc_matrix
from scipy.special import logsumexp, digamma, betaln
from .vireo_base import normalize, loglik_amplify, beta_entropy, get_binom_coeff


class BinomMixtureVB():
    """Binomial mixture model with variational inference

    The prior can be set via set_prior() before fitting the model.

    Key properties
    --------------
    beta_mu: numpy array (n_var, n_donor)
        Beta mean parameter of theta's posterior
    beta_sum: numpy array (n_var, n_donor)
        Beta concetration parameter of theta's posterior
    ID_prob: numpy array (n_cell, n_donor)
        Posterior cell assignment probability to each donor
    """

    def __init__(self, n_cell, n_var, n_donor, fix_beta_sum=False, 
        beta_mu_init=None, beta_sum_init=None, ID_prob_init=None):
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
        fix_beta_sum: bool. 
            Whether fixing the concetration parameter of theta's posterior
        beta_mu_init: numpy array (n_var, n_donor)
            Initial value of beta_mu, the mean parameter of theta
        beta_sum_init: numpy array (n_var, n_donor)
            Initial value of beta_sum, the concetration parameter of theta
        ID_prob_init: numpy array (n_cell, n_donor)
            Initial value of ID_prob, cell assignment probability to each donor
        """
        self.n_var = n_var
        self.n_cell = n_cell
        self.n_donor = n_donor
        self.fix_beta_sum = fix_beta_sum
        self.ID_prob_init = ID_prob_init
        self.beta_mu_init = beta_mu_init
        self.beta_sum_init = beta_sum_init
        
        # set priors; you can re-set by run this function
        self.set_prior()

        # initial key parameters
        self.set_initial(
            self.beta_mu_init, self.beta_sum_init, self.ID_prob_init
        )
        

    def set_initial(self, beta_mu_init=None, beta_sum_init=None, 
        ID_prob_init=None):
        """Random initialization
        """
        # initial key parameters
        if beta_mu_init is not None:
            self.beta_mu = beta_mu_init
        else:
            self.beta_mu = np.ones((self.n_var, self.n_donor)) * 0.5

        if beta_sum_init is not None:
            self.beta_sum = beta_sum_init
        else:
            self.beta_sum = np.ones(self.beta_mu.shape) * 30

        if ID_prob_init is not None:
            self.ID_prob = normalize(ID_prob_init, axis=1)
        else:
            self.ID_prob = normalize(np.random.rand(self.n_cell, self.n_donor))

        self.ELBO_iters = np.array([])
    
    def set_prior(self, ID_prior=None, beta_mu_prior=None, 
        beta_sum_prior=None):
        """Set prior for key variables: theta and ID_prob.
        The priors are in the same shape as its according variables.
        """
        if beta_mu_prior is None:
            beta_mu_prior = np.ones((self.n_var, self.n_donor)) * 0.5
        if beta_sum_prior is None:
            beta_sum_prior = np.ones(beta_mu_prior.shape) * 2.0

        self.theta_s1_prior = beta_mu_prior * beta_sum_prior
        self.theta_s2_prior = (1 - beta_mu_prior) * beta_sum_prior

        if ID_prior is not None:
            if len(ID_prior.shape) == 1:
                ID_prior = np.expand_dims(ID_prior, axis=0)
            self.ID_prior = ID_prior
        else:
            self.ID_prior = normalize(np.ones((self.n_cell, self.n_donor)))

    @property
    def theta_s1(self):
        """Beta concetration1 parameter for theta posterior"""
        return self.beta_mu * self.beta_sum

    @property
    def theta_s2(self):
        """Beta concetration2 parameter for theta posterior"""
        return (1 - self.beta_mu) * self.beta_sum


    def get_E_logLik(self, AD, DP):
        """Get the expecation of logLikelihood
        E_theta [P(AD|DP, theta, Z)]
        """
        BD = DP - AD

        # shape: (n_cell, n_donor)
        _E_logLik_mat = (
            AD.T @ digamma(self.theta_s1) + 
            BD.T @ digamma(self.theta_s2) -
            DP.T @ digamma(self.theta_s1 + self.theta_s2)
        )
        return _E_logLik_mat


    def update_theta_size(self, AD, DP):
        """Coordinate ascent for updating theta posterior parameters
        """
        BD = DP - AD
        _theta_s1 = AD @ self.ID_prob  #(n_var, n_donor)
        _theta_s2 = BD @ self.ID_prob  #(n_var, n_donor)
        _theta_s1 += self.theta_s1_prior
        _theta_s2 += self.theta_s2_prior

        self.beta_mu = _theta_s1 / (_theta_s1 + _theta_s2)
        if self.fix_beta_sum == False:
            self.beta_sum = _theta_s1 + _theta_s2


    def update_ID_prob(self, AD=None, DP=None, logLik_ID=None):
        """Coordinate ascent for updating assignment probability
        """
        if logLik_ID is None:
            logLik_ID = self.get_E_logLik(AD, DP)

        self.ID_prob = normalize(np.exp(loglik_amplify(
            logLik_ID + np.log(self.ID_prior))))                


    def get_ELBO(self, AD=None, DP=None, logLik_ID=None):
        """Calculating variational evidence lower bound with current parameters

        logLik_ID: numpy array (n_cell, n_donor), the output from update_ID_prob
        """
        if logLik_ID is None:
            self.get_E_logLik(AD, DP)

        LB_p = np.sum(logLik_ID * self.ID_prob)
        KL_ID = np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1))
        KL_theta = beta_entropy(
            np.append(
                np.expand_dims(self.theta_s1, 1), 
                np.expand_dims(self.theta_s2, 1), axis = 1),
            np.append(
                np.expand_dims(self.theta_s1_prior, 1), 
                np.expand_dims(self.theta_s2_prior, 1), axis = 1))
        
        return LB_p - KL_ID - KL_theta


    def _fit_BV(self, AD, DP, max_iter=200, min_iter=20, epsilon_conv=1e-2,
        verbose=True):
        """Fit Vireo model with coordinate ascent
        """
        ELBO = np.zeros(max_iter)
        for it in range(max_iter):
            self.update_theta_size(AD, DP)
            _logLik_ID = self.get_E_logLik(AD, DP)

            self.update_ID_prob(logLik_ID = _logLik_ID)
            ELBO[it] = self.get_ELBO(logLik_ID = _logLik_ID)

            if it > min_iter:
                if ELBO[it] - ELBO[it - 1] < -1e-6:
                    if verbose:
                        print("Warning: ELBO decreases %.8f to %.8f!\n"
                              %(ELBO[it - 1], ELBO[it]))
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif ELBO[it] - ELBO[it - 1] < epsilon_conv:
                    break

        self.ELBO_iters = np.append(self.ELBO_iters, ELBO[:it])


    def fit(self, AD, DP, n_init=10, max_iter=200, max_iter_pre=100, 
        random_seed=None, **kwargs):
        """Fit VB with multiple initializations

        Parameters
        ----------
        AD : scipy.sparse.csc_matrix (n_var, n_cell)
            Sparse count matrix for alternative allele
        DP : scipy.sparse.csc_matrix (n_var, n_cell)
            Sparse count matrix for depths, alternative + refeerence alleles
        n_inits : int
            Number of random initialisations to use
        max_iter : int
            Maximum number of iterations for _fit_BV() in best initial
        max_iter_pre : int
            Maximum number of iterations for _fit_BV() in multiple initials
        min_iter :
            Minimum number of iterations for _fit_BV()
        epsilon_conv : float
            Threshold for detecting convergence for _fit_BV()
        verbose : bool
            Whether print out log info for _fit_BV()
        random_seed : None or int
            Random seed in numpy.random for multiple initializations
        """
        if random_seed is not None:
            np.random.seed(random_seed)
            
        if type(DP) is np.ndarray and np.mean(DP > 0) < 0.3:
            print("Warning: input matrices is %.1f%% sparse, " 
                  %(100 - np.mean(DP > 0) * 100) +
                  "change to scipy.sparse.csc_matrix" )
            AD = csc_matrix(AD)
            DP = csc_matrix(DP)

        _binom_coeff = np.sum(get_binom_coeff(AD, DP, is_log = True))

        self.ELBO_inits = []
        for i in range(n_init):
            self.set_initial(
                self.beta_mu_init, self.beta_sum_init, self.ID_prob_init
            )
            self._fit_BV(AD, DP, max_iter=max_iter_pre, **kwargs)
            self.ELBO_inits.append(self.ELBO_iters[-1])
            
            ## first or better initialization
            if i == 0 or (self.ELBO_iters[-1] > np.max(self.ELBO_inits[:-1])):
                _ID_prob_best = self.ID_prob + 0
                _beta_mu_best = self.beta_mu + 0
                _beta_sum_best = self.beta_sum + 0
                _ELBO_iters_best = self.ELBO_iters + 0
                
        ## Re-fit with best parameters
        self.set_initial(_beta_mu_best, _beta_sum_best, _ID_prob_best)
        self.ELBO_iters = _ELBO_iters_best
        self._fit_BV(AD, DP, max_iter=max_iter, **kwargs)

        ## add binomial coefficient constants
        self.ELBO_iters = self.ELBO_iters + _binom_coeff
        self.ELBO_inits = np.array(self.ELBO_inits) + _binom_coeff
        