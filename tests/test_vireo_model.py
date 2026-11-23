import unittest
from math import comb, log

import numpy as np
from scipy.sparse import csc_matrix
from scipy.special import digamma, softmax

from vireoSNP import Vireo


class VireoModelTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(123)
        allele_fractions = np.array([[0.02, 0.98], [0.5, 0.98], [0.98, 0.02], [0.98, 0.5]])
        donors = np.repeat([0, 1], 8)
        self.depth = rng.integers(4, 18, (4, donors.size))
        self.depth[::2, ::3] = 0
        self.counts = rng.binomial(self.depth, allele_fractions[:, donors])
        self.ad = csc_matrix(self.counts)
        self.dp = csc_matrix(self.depth)
        self.prior = np.array([[0.5, 0.3, 0.2], [0.1, 0.3, 0.6]])

    def make_model(self, prior=None, **kwargs):
        n_var, n_cell = self.depth.shape
        n_donor, n_gt = self.prior.shape
        rng = np.random.default_rng(71)
        model = Vireo(
            n_cell=n_cell,
            n_var=n_var,
            n_donor=n_donor,
            n_GT=n_gt,
            ID_prob_init=rng.random((n_cell, n_donor)),
            GT_prob_init=rng.random((n_var, n_donor, n_gt)),
            **kwargs
        )
        if prior is not None:
            model.set_prior(GT_prior=prior)
        return model

    def fit_at_limit(self, model, iterations, sparse=True):
        ad, dp = (self.ad, self.dp) if sparse else (self.counts, self.depth)
        model.fit(ad, dp, max_iter=iterations, min_iter=iterations, verbose=False)

    def assert_same_fit(self, first, second):
        for attr in ("GT_prob", "ID_prob", "beta_mu", "beta_sum", "ELBO_"):
            np.testing.assert_allclose(
                getattr(first, attr), getattr(second, attr), rtol=0, atol=1e-12
            )

    def assert_final_bound_matches_posteriors(self, model):
        # fit() adds the log-binomial coefficient once to each computed bound.
        terms = [log(comb(int(dp), int(ad))) for ad, dp in zip(self.counts.flat, self.depth.flat)]
        expected = model.get_ELBO(None, AD=self.ad, DP=self.dp) + sum(terms)
        # get_binom_coeff stores this data-only term in float32.
        np.testing.assert_allclose(model.ELBO_[-1], expected, rtol=0, atol=1e-4)

    def test_csc_fit_with_broadcast_prior_matches_expanded_prior(self):
        prior_before_fit = self.prior.copy()
        expanded = np.broadcast_to(self.prior, (4, 2, 3)).copy()
        full_model = self.make_model(expanded)
        broadcast_model = self.make_model(self.prior)
        self.fit_at_limit(full_model, 10)
        self.fit_at_limit(broadcast_model, 10)

        self.assert_same_fit(broadcast_model, full_model)
        self.assertEqual(broadcast_model.GT_prior.shape, (1, 2, 3))
        self.assertEqual(broadcast_model.GT_prob.shape, (4, 2, 3))
        self.assertEqual(broadcast_model.ID_prob.shape, (16, 2))
        np.testing.assert_allclose(broadcast_model.GT_prob.sum(axis=-1), 1)
        np.testing.assert_allclose(broadcast_model.ID_prob.sum(axis=-1), 1)
        np.testing.assert_array_equal(broadcast_model.ID_prob.argmax(axis=1),
                                      full_model.ID_prob.argmax(axis=1))
        np.testing.assert_array_equal(self.prior, prior_before_fit)

    def test_sparse_genotype_step_matches_scalar_posterior_oracle(self):
        model = self.make_model(self.prior, ASE_mode=True)
        expected = np.empty_like(model.GT_prob)
        for variant, (ad, dp) in enumerate(zip(self.counts, self.depth)):
            alpha = model.beta_mu[variant] * model.beta_sum[variant]
            beta = (1 - model.beta_mu[variant]) * model.beta_sum[variant]
            alt_log = digamma(alpha) - digamma(alpha + beta)
            ref_log = digamma(beta) - digamma(alpha + beta)
            for donor, prior in enumerate(self.prior):
                scores = np.log(prior).copy()
                for cell, (alternate, depth) in enumerate(zip(ad, dp)):
                    scores += model.ID_prob[cell, donor] * (
                        alternate * alt_log + (depth - alternate) * ref_log
                    )
                expected[variant, donor] = softmax(scores)

        model.update_GT_prob(self.ad, self.dp)
        np.testing.assert_allclose(model.GT_prob, expected, rtol=0, atol=1e-12)

    def test_variant_specific_full_prior_agrees_on_dense_and_sparse_counts(self):
        full_prior = np.stack([np.roll(self.prior, i, axis=1) for i in range(4)])
        sparse_model = self.make_model(full_prior.copy())
        dense_model = self.make_model(full_prior.copy())
        self.fit_at_limit(sparse_model, 6)
        self.fit_at_limit(dense_model, 6, sparse=False)
        self.assert_same_fit(sparse_model, dense_model)
        np.testing.assert_allclose(sparse_model.GT_prior, full_prior, rtol=0, atol=1e-12)

    def test_default_prior_matches_explicit_uniform_prior(self):
        default_model = self.make_model()
        explicit_model = self.make_model(np.full((4, 2, 3), 1 / 3))
        self.fit_at_limit(default_model, 6)
        self.fit_at_limit(explicit_model, 6)
        self.assert_same_fit(default_model, explicit_model)

    def test_iteration_limit_preserves_seeded_history_prefix_and_final_bound(self):
        ten, eleven = self.make_model(), self.make_model()
        self.fit_at_limit(ten, 10)
        self.fit_at_limit(eleven, 11)
        self.assertEqual(ten.ELBO_.size, 10)
        self.assertEqual(eleven.ELBO_.size, 11)
        np.testing.assert_allclose(ten.ELBO_, eleven.ELBO_[:10], rtol=0, atol=1e-12)
        self.assertTrue(np.isfinite(eleven.ELBO_).all())
        self.assert_final_bound_matches_posteriors(ten)
        self.assert_final_bound_matches_posteriors(eleven)

    def test_early_convergence_keeps_the_converging_step(self):
        model = self.make_model(learn_GT=False, learn_theta=False)
        # With both held fixed, the assignment update reaches a fixed point.
        model.fit(self.ad, self.dp, max_iter=8, min_iter=0, verbose=False)
        self.assertEqual(model.ELBO_.size, 2)
        np.testing.assert_allclose(model.ELBO_[0], model.ELBO_[1], rtol=0, atol=1e-12)
        self.assert_final_bound_matches_posteriors(model)

    def test_one_iteration_keeps_its_only_bound(self):
        model = self.make_model()
        self.fit_at_limit(model, 1)
        self.assertEqual(model.ELBO_.size, 1)
        self.assertTrue(np.isfinite(model.ELBO_[0]))
        self.assert_final_bound_matches_posteriors(model)


if __name__ == "__main__":
    unittest.main()
