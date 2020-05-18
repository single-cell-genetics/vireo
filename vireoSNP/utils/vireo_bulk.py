# Identification of donor abundance in bulk sample
import numpy as np

__docformat__ = "restructuredtext en"

__all__ = ['VireoBulk']

class VireoBulk():
    """
    Estimate of donor abundance in a multipexed bulk sample

    Varibale to infer
    -----------------
    psi: numpy.array (n_donor, )
        The fractional abundance of each donor in the mixture
    theta: numpy.array (n_GT, )
        The alternative allele rate in each genotype category

    Parameters
    ----------
    n_GT: int, number of genotype categories
    n_donor: int, number of donors in the mixture
    """
    def __init__(self, n_donor, n_GT=3, psi_init=None, 
                 theta_init=[0.01, 0.5, 0.99]):
        self.n_GT = n_GT
        self.n_donor = n_donor
        
        self.psi = np.random.dirichlet([1] * n_donor)
        self.theta = np.random.rand(n_GT)
        
        if psi_init is not None:
            if n_donor != len(psi_init):
                print("Warning: n_donor != len(psi_init)")
            else:
                self.psi = np.random.dirichlet([1] * n_donor)

        if theta_init is not None:
            if n_GT != len(theta_init):
                print("Warning: n_GT != len(theta_init)")
            else:
                self.theta = theta_init
    
    def fit(self, AD, DP, GT_prob, max_iter=200, min_iter=5, epsilon_conv=1e-3,
            learn_theta=True, delay_fit_theta=0, model="EM", verbose=False):
        """Fit the unknown variable psi and theta with EM algorithm

        Parameters
        ----------
        AD: numpy.array, (n_variant, ), int
            The count vector for alternative allele in all variants
        DP: numpy.array (n_variant, ), int
            The count vector for depths in all variants (i.e., two alleles)
        GT_prob: numpy.array, (n_variants, n_donor, n_GT)
            The probability tensor for each genotype in each donor
        learn_theta: bool
            Whether learn theta, otherwise use theta_init
        delay_fit_theta: int
            The number of steps to delay in updating theta
        max_iter : int
            Maximum number of iterations
        min_iter :
            Minimum number of iterations
        epsilon_conv : float
            Threshold for detecting convergence
        model: string
            The algorithm used to fit the model. Only "EM" is supported for
            Expectation-Maximumization algorithm
        verbose : bool
            Whether print out log info
        """
        BD = DP - AD
        
        logLik = np.zeros(max_iter)
        for it in range(max_iter):
            # E step: expectation of count assignment probability
            theta_mat = np.tensordot(GT_prob, self.theta, axes=(2, 0)) #np.dot
            Z1 = theta_mat * np.expand_dims(self.psi, 0)
            Z1 = Z1 / np.sum(Z1, axis=1, keepdims=True)
            
            Z0 = (1 - theta_mat) * np.expand_dims(self.psi, 0)
            Z0 = Z0 / np.sum(Z0, axis=1, keepdims=True)

            # M step: maximize logLikehood over psi and theta
            psi_raw = np.dot(AD, Z1) + np.dot(BD, Z0)
            self.psi = psi_raw / np.sum(psi_raw)
            
            if learn_theta and it >= delay_fit_theta:
                theta_s1 = np.dot(AD, np.sum(GT_prob * np.expand_dims(Z1, 2), axis = 1))
                theta_s2 = np.dot(BD, np.sum(GT_prob * np.expand_dims(Z0, 2), axis = 1))
                self.theta = theta_s1 / (theta_s1 + theta_s2)
            
            # Likelihood and check convergence
            theta_vct = np.dot(np.dot(GT_prob, self.theta), self.psi)
            logLik[it] = np.sum(
                AD * np.log(theta_vct) + BD * np.log(1 - theta_vct))
            if it > min_iter:
                if logLik[it] < logLik[it - 1]:
                    if verbose:
                        print("Warning: logLikelihood decreases!\n")
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif logLik[it] - logLik[it - 1] < epsilon_conv:
                    break
        
        self.logLik = logLik[it]
        self.logLik_all = logLik[:it]

    def LR_test(self, **kwargs):
        """Likelihood ratio test for psi vector in a null hypothesis.

        Use **kwargs for psi_null, AD, DP, GT_prob, log in 
        vireoSNP.LikRatio_test() function. Note, AD, DP, GT_prob the same as 
        the self.fit() function.
        """
        return LikRatio_test(psi=self.psi, theta=self.theta, **kwargs)


def LikRatio_test(psi, psi_null, AD, DP, GT_prob, theta, log=False):
    """Likelihood ratio test for psi vector in a null hypothesis.
    Please use the same AD, DP, and GT_prob as the fit() function.

    Parameters
    ----------
    psi: numpy.array (n_donor, )
        The fractional abundance of each donor in the mixture for alternative
        hypothesis
    psi_null: numpy.array (n_donor, )
        The psi vector in a null hypothesis
    AD: numpy.array, (n_variant, ), int
        The count vector for alternative allele in all variants
    DP: numpy.array (n_variant, ), int
        The count vector for depths in all variants (i.e., two alleles)
    GT_prob: numpy.array, (n_variants, n_donor, n_GT)
        The probability tensor for each genotype in each donor
    theta: numpy.array (n_GT, )
        The alternative allele rate in each genotype category
    log: bool
        If return p value in logarithm scale

    Return
    ------
    statistic: float
        The calculated chi2-statistic.
    pvalue: float
        The single-tailed p-value.
    """
    from scipy.stats import chi2

    BD = DP - AD
    theta_vct_alt = np.dot(np.dot(GT_prob, theta), psi)
    logLik_alt = np.sum(
        AD * np.log(theta_vct_alt) + BD * np.log(1 - theta_vct_alt))

    theta_vct_null = np.dot(np.dot(GT_prob, theta), psi_null)
    logLik_null = np.sum(
        AD * np.log(theta_vct_null) + BD * np.log(1 - theta_vct_null))

    LR = 2 * (logLik_alt - logLik_null)
    df = len(psi_null) - 1
    if log:
        pval = chi2.logsf(LR, df)
    else:
        pval = chi2.sf(LR, df)
    
    return LR, pval
