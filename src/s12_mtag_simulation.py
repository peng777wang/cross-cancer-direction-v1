"""Simulation: does MTAG manufacture opposite-direction pleiotropic loci?

Two traits, true effects drawn from a bivariate normal with genetic correlation
rho; sample sizes N_A (large) and N_B (small).  MTAG is applied with Omega
estimated from the simulated data, exactly as in practice.

Result (see docs/METHODS_NOTES.md): under a correctly specified model MTAG never
flips a sign -- every opposite-direction locus it reports is genuinely opposite.
What it does do is *amplify* a weak trait's estimate by borrowing from the strong
trait, by up to ~13% at rho=-0.5 with a 20-fold power difference.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from s10_mtag_core import estimate_omega_mom, mtag_estimate

SIG = 5.451


def run(rho, ratio, seed, M, N_A, h):
    rng = np.random.default_rng(seed)
    N = np.array([N_A, N_A * ratio], dtype=float)
    cov = h * np.array([[1.0, rho], [rho, 1.0]])
    b = rng.multivariate_normal([0, 0], cov, size=M)
    Z = np.sqrt(N)[None, :] * b + rng.standard_normal((M, 2))
    om = estimate_omega_mom(Z, np.tile(N, (M, 1)))
    om = (om + om.T) / 2
    mt, se = mtag_estimate(Z, np.tile(N, (M, 1)), om, np.eye(2))
    mtz = mt / se
    both = (np.abs(mtz) > SIG).all(axis=1)
    opp = both & (np.sign(mtz[:, 0]) != np.sign(mtz[:, 1]))
    false_opp = opp & (np.sign(b[:, 0]) == np.sign(b[:, 1]))
    with np.errstate(invalid="ignore", divide="ignore"):
        borrow = 1 - np.abs(Z[both, 1]) / np.abs(mtz[both, 1])
    return {"rho": rho, "N_ratio": ratio, "n_both_sig": int(both.sum()),
            "n_mtag_opposite": int(opp.sum()),
            "n_opposite_but_true_concordant": int(false_opp.sum()),
            "median_borrow_traitB": float(np.nanmedian(borrow)) if both.sum() else np.nan}


def main(seed0, reps, M, N_A, h, out):
    rhos = [-0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5]
    ratios = [1.0, 0.2, 0.05]
    rows = []
    for i_rho, rho in enumerate(rhos):
        for i_ratio, ratio in enumerate(ratios):
            for rep in range(reps):
                rows.append(run(rho, ratio, seed0 + 1000 * rep + 37 * i_rho + i_ratio,
                                M, N_A, h))
    os.makedirs(out, exist_ok=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, "mtag_simulation_raw.csv"), index=False)
    g = res.groupby(["rho", "N_ratio"]).agg(
        both_sig=("n_both_sig", "mean"),
        opposite=("n_mtag_opposite", "mean"),
        false_opposite=("n_opposite_but_true_concordant", "mean"),
        median_borrow=("median_borrow_traitB", "mean")).reset_index()
    g.to_csv(os.path.join(out, "mtag_simulation_summary.csv"), index=False)
    print(g.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--snps", type=int, default=200000)
    ap.add_argument("--n-large", type=int, default=100000)
    ap.add_argument("--h2-snp", type=float, default=2.0e-4)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    main(a.seed, a.reps, a.snps, a.n_large, a.h2_snp, a.out)
