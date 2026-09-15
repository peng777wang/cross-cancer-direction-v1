"""Exact port of coloc::combine.abf / coloc.abf (coloc v5.2.3).

Verified against the official R source (R/claudia.R, lines 109-133) and on four
synthetic scenarios (shared variant -> H4; one-sided -> H1; two distinct -> H3;
no signal -> H0).  The H3 prior is p1*p2, NOT p12 -- using p12 inflates H3 by
10,000x at the default priors and produced a false-negative screen early in this
project (docs/AUDIT_LOG.md, item 1).

    lABF_i = 0.5 * (log(1 - r) + r * z_i^2),   r = V_prior / (V_prior + se_i^2)
    H0 ~ 1
    H1 ~ p1  * sum exp(lABF1)
    H2 ~ p2  * sum exp(lABF2)
    H3 ~ p1*p2 * (sum e1 * sum e2 - sum e1*e2)
    H4 ~ p12 * sum e1*e2
"""
import numpy as np

from common import logdiff, logsumexp


def labf(beta, se, sd_prior):
    """Log approximate Bayes factor for one SNP (coloc::approx.bf.estimates)."""
    v = np.asarray(se, dtype=float) ** 2
    V = float(sd_prior) ** 2
    r = V / (V + v)
    z = np.asarray(beta, dtype=float) / np.asarray(se, dtype=float)
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def coloc_abf(beta1, se1, beta2, se2, sd1=0.15, sd2=0.2,
              p1=1e-4, p2=1e-4, p12=1e-5):
    """Posterior probabilities for H0..H4.  sd1/sd2 are the coloc defaults
    (0.15 for a quantitative trait, 0.2 for a case/control trait)."""
    l1, l2 = labf(beta1, se1, sd1), labf(beta2, se2, sd2)
    ls1, ls2 = logsumexp(l1), logsumexp(l2)
    ls12 = logsumexp(l1 + l2)
    lH = np.array([0.0,
                   np.log(p1) + ls1,
                   np.log(p2) + ls2,
                   np.log(p1) + np.log(p2) + logdiff(ls1 + ls2, ls12),
                   np.log(p12) + ls12])
    m = lH.max()
    pp = np.exp(lH - m)
    return pp / pp.sum(), int(len(l1))


if __name__ == "__main__":
    # the four verification scenarios
    rng = np.random.default_rng(0)
    n, se = 2000, np.full(2000, 0.08)
    b1 = np.zeros(n); b1[100] = 0.35
    b3 = np.zeros(n); b3[700] = 0.35
    zero = np.zeros(n)
    cases = {"same shared variant": (b1, b1), "trait 1 only": (b1, zero),
             "trait 2 only": (zero, b1), "two distinct variants": (b1, b3),
             "shared but flipped": (b1, -b1), "neither": (zero, zero)}
    for label, (x, y) in cases.items():
        pp, _ = coloc_abf(x, se, y, se)
        print("%-22s H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f" % ((label,) + tuple(pp)))
