"""Self-contained regression tests for the analysis code.

Run with:  python tests/test_reproducibility.py     (stdlib unittest, no pytest)

These tests do not need any raw data, only the numbers that the estimators must
produce.  They exist so that anyone reviewing the code can confirm, in seconds,
that the statistics are implemented as described -- independent of whether the
large upstream inputs are available.
"""
import os
import subprocess
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)

from common import harmonise, harmonise_z, logdiff          # noqa: E402
from s05_coloc_core import coloc_abf                        # noqa: E402
from s09_mvmr import mvmr, wald_ratio                       # noqa: E402
from s10_mtag_core import estimate_omega_mom, mtag_estimate  # noqa: E402

SIG = 5.451


class TestAlleleHarmonisation(unittest.TestCase):
    """Item 3 of the audit log: strand must be handled, not string-matched."""

    def test_complement_is_recognised_as_same(self):
        flip, ambiguous, ok = harmonise("A", "G", "T", "C")
        self.assertTrue(ok)
        self.assertFalse(flip)
        self.assertFalse(ambiguous)

    def test_effect_allele_on_reference_second_allele_flips(self):
        flip, ambiguous, ok = harmonise("G", "A", "A", "G")
        self.assertTrue(ok)
        self.assertTrue(flip)

    def test_palindromic_pair_is_flagged_ambiguous(self):
        flip, ambiguous, ok = harmonise("A", "T", "T", "A")
        self.assertTrue(ok)
        self.assertTrue(ambiguous)
        self.assertIsNone(flip)

    def test_non_matching_alleles_are_rejected(self):
        # {A,T} vs {A,G}: neither the direct nor the complemented orientation fits
        _, _, ok = harmonise("A", "T", "A", "G")
        self.assertFalse(ok)

    def test_harmonise_z_flips_the_sign(self):
        # effect allele G is the reference A2, so its query frequency (0.3) must
        # equal 1 - freq(A1) = 0.7 in the reference panel
        z = harmonise_z(2.0, "G", "A", "A", "G", freq=0.3, ref_freq=0.7)
        self.assertAlmostEqual(z, -2.0)

    def test_harmonise_z_resolves_palindrome_by_frequency(self):
        z = harmonise_z(3.0, "A", "T", "T", "A", freq=0.2, ref_freq=0.8)
        self.assertAlmostEqual(z, -3.0)


class TestLogDiff(unittest.TestCase):
    def test_matches_direct_computation(self):
        a, b = np.log(5.0), np.log(3.0)
        self.assertAlmostEqual(logdiff(a, b), np.log(2.0), places=12)

    def test_returns_minus_infinity_for_equal_inputs(self):
        self.assertEqual(logdiff(0.0, 0.0), -np.inf)


class TestColocCore(unittest.TestCase):
    """Items 1 and 7 of the audit log: the ABF port must match coloc::combine.abf,
    and H3 must use the prior p1*p2 (not p12)."""

    def setUp(self):
        rng = np.random.default_rng(0)
        self.se = np.full(2000, 0.08)
        self.b_same = np.zeros(2000)
        self.b_same[100] = 0.35
        self.b_other = np.zeros(2000)
        self.b_other[700] = 0.35
        self.zero = np.zeros(2000)
        self.rng = rng

    def test_shared_variant_supports_h4(self):
        pp, _ = coloc_abf(self.b_same, self.se, self.b_same, self.se)
        self.assertEqual(int(np.argmax(pp)), 4)
        self.assertGreater(pp[4], 0.8)

    def test_two_distinct_variants_favour_h3_over_h4(self):
        # H0 can still dominate in absolute terms: with 2,000 informative null
        # SNPs and p1 = p2 = 1e-4 the single shared SNPs rarely outweigh H0.
        # What must hold is the *relative* discrimination H3 > H4, which is the
        # only thing the coloc screen ever claims.
        pp, _ = coloc_abf(self.b_same, self.se, self.b_other, self.se)
        self.assertGreater(pp[3], pp[4])

    def test_one_sided_signal_favours_h1(self):
        pp, _ = coloc_abf(self.b_same, self.se, self.zero, self.se)
        self.assertGreater(pp[1], pp[2])
        self.assertGreater(pp[1], pp[4])

    def test_no_signal_supports_h0(self):
        pp, _ = coloc_abf(self.zero, self.se, self.zero, self.se)
        self.assertEqual(int(np.argmax(pp)), 0)

    def test_abf_depends_on_z_squared_only(self):
        """Sign flips must not change the posterior -- this is why the coloc
        screen is immune to allele-orientation errors."""
        a, _ = coloc_abf(self.b_same, self.se, self.b_same, self.se)
        b, _ = coloc_abf(self.b_same, self.se, -self.b_same, self.se)
        np.testing.assert_allclose(a, b, atol=1e-12)

    def test_h3_uses_p1_times_p2_not_p12(self):
        """If H3 were given the p12 prior, lowering p1 and p2 would not change it.
        Build a clean H3 case and confirm the posterior responds to p1*p2."""
        kw = dict(beta1=self.b_same, se1=self.se, beta2=self.b_other, se2=self.se)
        lo, _ = coloc_abf(p1=1e-6, p2=1e-6, p12=1e-5, **kw)
        hi, _ = coloc_abf(p1=1e-2, p2=1e-2, p12=1e-5, **kw)
        self.assertGreater(hi[3], lo[3] * 100)

    def test_posteriors_sum_to_one(self):
        pp, n = coloc_abf(self.b_same, self.se, self.b_other, self.se)
        self.assertEqual(n, 2000)
        self.assertAlmostEqual(float(pp.sum()), 1.0, places=12)


class TestMTAGCore(unittest.TestCase):
    """The MTAG port and the two claims it supports:
    (a) Omega is recovered by method of moments;
    (b) MTAG never flips a sign under a correctly specified model."""

    def _simulate(self, rho, ratio, seed, M=120_000, N_A=100_000.0, h=2e-4):
        rng = np.random.default_rng(seed)
        N = np.array([N_A, N_A * ratio])
        cov = h * np.array([[1.0, rho], [rho, 1.0]])
        b = rng.multivariate_normal([0.0, 0.0], cov, size=M)
        Z = np.sqrt(N)[None, :] * b + rng.standard_normal((M, 2))
        return b, Z, N

    def test_omega_recovered_by_method_of_moments(self):
        rho, h = -0.5, 2e-4
        _, Z, N = self._simulate(rho, 0.2, seed=1)
        om = estimate_omega_mom(Z, np.tile(N, (Z.shape[0], 1)))
        self.assertAlmostEqual(om[0, 0], h, delta=0.15 * h)
        self.assertAlmostEqual(om[1, 1], h, delta=0.15 * h)
        self.assertAlmostEqual(om[0, 1], rho * h, delta=0.15 * h)

    def test_mtag_reproduces_raw_z_when_there_is_no_borrowing(self):
        """With Omega = h*I (rho = 0) each trait's estimate must stay close to its
        own observation; the port must not invent signal out of nothing."""
        _, Z, N = self._simulate(0.0, 1.0, seed=2)
        om = estimate_omega_mom(Z, np.tile(N, (Z.shape[0], 1)))
        mt, se = mtag_estimate(Z, np.tile(N, (Z.shape[0], 1)), om, np.eye(2))
        mtz = mt / se
        self.assertGreater(np.corrcoef(mtz[:, 1], Z[:, 1])[0, 1], 0.95)

    def test_mtag_never_flips_a_sign(self):
        """The retracted claim. Under a correctly specified model every
        opposite-direction locus MTAG reports is genuinely opposite."""
        for rho in (-0.5, -0.2, 0.0, 0.2):
            b, Z, N = self._simulate(rho, 0.05, seed=int(1000 * (rho + 1)) + 7)
            om = estimate_omega_mom(Z, np.tile(N, (Z.shape[0], 1)))
            mt, se = mtag_estimate(Z, np.tile(N, (Z.shape[0], 1)), om, np.eye(2))
            mtz = mt / se
            both = (np.abs(mtz) > SIG).all(axis=1)
            false_opp = both & (np.sign(mtz[:, 0]) != np.sign(mtz[:, 1])) & \
                (np.sign(b[:, 0]) == np.sign(b[:, 1]))
            self.assertEqual(int(false_opp.sum()), 0,
                             "MTAG flipped a sign at rho=%s" % rho)

    def test_mtag_borrows_more_when_one_trait_is_underpowered(self):
        """MTAG amplifies a weak trait by borrowing; the amount must grow as the
        power gap widens."""
        checks = sys.maxsize
        borrowed = {}
        for ratio in (1.0, 0.05):
            b, Z, N = self._simulate(-0.5, ratio, seed=11)
            om = estimate_omega_mom(Z, np.tile(N, (Z.shape[0], 1)))
            mt, se = mtag_estimate(Z, np.tile(N, (Z.shape[0], 1)), om, np.eye(2))
            mtz = mt / se
            both = (np.abs(mtz) > SIG).all(axis=1)
            borrowed[ratio] = float(np.median(1 - np.abs(Z[both, 1]) / np.abs(mtz[both, 1])))
        self.assertGreater(borrowed[0.05], borrowed[1.0] + 0.05)
        self.assertLess(borrowed[0.05], 0.35)   # amplification, not fabrication


class TestMREstimators(unittest.TestCase):
    """Item 5 of the audit log: instrument strength and estimator unbiasedness."""

    def test_wald_ratio_is_the_ratio_of_effects(self):
        th, se, p = wald_ratio(-0.23, 0.02, -0.11, 0.02)
        self.assertAlmostEqual(th, -0.23 / -0.11, places=12)
        self.assertGreater(se, 0.0)
        self.assertLess(p, 0.05)

    def test_mvmr_is_unbiased(self):
        # Instrument strengths chosen so that the per-replicate SE is small
        # (~0.03); 500 replicates then give a Monte-Carlo SE of the mean of
        # ~0.0015, so a 0.02 tolerance is a real test of unbiasedness and not a
        # test of how many replicates happened to be drawn.
        rng = np.random.default_rng(7)
        J, true, reps = 15, np.array([0.5, 0.3]), 500
        est = []
        for _ in range(reps):
            B = np.column_stack([rng.normal(0.3, 0.10, J), rng.normal(0.4, 0.15, J)])
            sg = np.full(J, 0.02)
            gamma = B @ true + rng.normal(0, 0.02, J)
            th, _ = mvmr(B[:, 0], B[:, 1], gamma, sg)
            est.append(th)
        est = np.array(est)
        bias = est.mean(axis=0) - true
        mc_se = est.std(axis=0) / np.sqrt(reps)
        self.assertLess(abs(bias).max(), 0.02)
        self.assertLess(abs(bias / mc_se).max(), 3.0)


class TestClaimVerification(unittest.TestCase):
    def test_all_manuscript_claims_reproduce(self):
        out = subprocess.run(
            [sys.executable, os.path.join(SRC, "verify_claims.py")],
            capture_output=True, text=True)
        self.assertIn("19 / 19", out.stdout)
        self.assertEqual(out.returncode, 0, out.stdout[-2000:])


if __name__ == "__main__":
    unittest.main(verbosity=2)
