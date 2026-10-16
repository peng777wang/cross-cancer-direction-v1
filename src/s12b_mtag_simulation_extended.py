"""Extended simulation: does the sign-safety of MTAG hold in the regime the
manuscript actually analyses?

Two additions over s12_mtag_simulation.py.
1. Sample-size ratios down to 1:200, covering the observed 1:77 disparity.
2. A shared-control scenario, in which the sampling errors of the two traits
   are correlated (rho_e > 0) while Omega is still estimated as if they were
   independent, which is what happens when case-control studies share controls.

Reported per scenario: number of loci significant in both traits, number that
MTAG calls opposite-direction, how many of those were truly concordant in the
generating model (these are the false opposites), and the median borrowing
fraction of the weak trait.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE, os.path.join(HERE, "..", "release", "src")):
    if os.path.exists(os.path.join(cand, "s10_mtag_core.py")):
        sys.path.insert(0, os.path.abspath(cand))
        break
from s10_mtag_core import estimate_omega_mom, mtag_estimate

SIG = 5.451


def run(rho, ratio, rho_e, seed, M, N_A, h):
    rng = np.random.default_rng(seed)
    N = np.array([N_A, N_A * ratio], dtype=float)
    cov = h * np.array([[1.0, rho], [rho, 1.0]])
    b = rng.multivariate_normal([0, 0], cov, size=M)
    err_cov = np.array([[1.0, rho_e], [rho_e, 1.0]])
    e = rng.multivariate_normal([0, 0], err_cov, size=M)
    Z = np.sqrt(N)[None, :] * b + e
    om = estimate_omega_mom(Z, np.tile(N, (M, 1)))
    om = (om + om.T) / 2
    mt, se = mtag_estimate(Z, np.tile(N, (M, 1)), om, np.eye(2))
    mtz = mt / se
    both = (np.abs(mtz) > SIG).all(axis=1)
    opp = both & (np.sign(mtz[:, 0]) != np.sign(mtz[:, 1]))
    false_opp = opp & (np.sign(b[:, 0]) == np.sign(b[:, 1]))
    with np.errstate(invalid="ignore", divide="ignore"):
        borrow = 1 - np.abs(Z[both, 1]) / np.abs(mtz[both, 1])
    return {"rho": rho, "N_ratio": ratio, "rho_e": rho_e,
            "n_both_sig": int(both.sum()),
            "n_mtag_opposite": int(opp.sum()),
            "n_false_opposite": int(false_opp.sum()),
            "median_borrow_traitB": float(np.nanmedian(borrow)) if both.sum()
            else np.nan}


def main(seed0, reps, M, N_A, h, out):
    rows = []
    ratios = [1.0, 0.2, 0.05, 0.02, 0.013, 0.005]
    rhos = [-0.5, -0.25, 0.0, 0.25, 0.5]
    for i_ratio, ratio in enumerate(ratios):
        for i_rho, rho in enumerate(rhos):
            for rep in range(reps):
                rows.append(run(rho, ratio, 0.0,
                                seed0 + 1000 * rep + 37 * i_rho
                                + 11 * i_ratio, M, N_A, h))
    cat = {"scenario": "independent traits"}
    for rho_e in (0.1, 0.3, 0.5):
        for i_ratio, ratio in enumerate((0.05, 0.013)):
            for i_rho, rho in enumerate((-0.5, 0.0)):
                for rep in range(reps):
                    rows.append(run(rho, ratio, rho_e,
                                    seed0 + 5000 + 1000 * rep + 37 * i_rho
                                    + 11 * i_ratio + int(100 * rho_e),
                                    M, N_A, h))
    res = pd.DataFrame(rows)
    os.makedirs(out, exist_ok=True)
    res.to_csv(os.path.join(out, "mtag_simulation_extended_raw.csv"),
               index=False)
    g = res.groupby(["rho_e", "rho", "N_ratio"]).agg(
        both_sig=("n_both_sig", "mean"),
        opposite=("n_mtag_opposite", "mean"),
        false_opposite=("n_false_opposite", "mean"),
        median_borrow=("median_borrow_traitB", "mean")).reset_index()
    g.to_csv(os.path.join(out, "mtag_simulation_extended_summary.csv"),
             index=False)
    pd.set_option("display.width", 200)
    print(g.round(3).to_string(index=False))
    return g


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261013)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--snps", type=int, default=200000)
    ap.add_argument("--n-large", type=int, default=100000)
    ap.add_argument("--h2-snp", type=float, default=2.0e-4)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    main(a.seed, a.reps, a.snps, a.n_large, a.h2_snp, a.out)
