"""Validate the MTAG port against a real MTAG run.

Given a chunk of MTAG output (Z, N, mtag_z) and the run's omega/sigma matrices,
this reproduces mtag_z from the port.  Expected: Pearson r = 1.00000000 and
maximum absolute difference ~1e-15.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from s10_mtag_core import mtag_estimate


def main(omega_txt, sigma_txt, mtag_dir, traits, n_snps, out_dir):
    om = np.loadtxt(omega_txt)
    sg = np.loadtxt(sigma_txt)
    frames = {}
    for i, t in enumerate(traits, start=1):
        d = pd.read_csv(os.path.join(mtag_dir, "mtag_all_trait_%d.txt" % i), sep="\t",
                        usecols=["SNP", "Z", "N", "mtag_z"], low_memory=False)
        frames[t] = d.iloc[:n_snps].set_index("SNP")
    snp = frames[traits[0]].index
    Z = np.column_stack([frames[t].loc[snp, "Z"].values for t in traits])
    N = np.column_stack([frames[t].loc[snp, "N"].values for t in traits])
    ref = np.column_stack([frames[t].loc[snp, "mtag_z"].values for t in traits])
    b, s = mtag_estimate(Z, N, om, sg)
    mine = b / s
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for j, t in enumerate(traits):
        r = np.corrcoef(mine[:, j], ref[:, j])[0, 1]
        md = float(np.nanmax(np.abs(mine[:, j] - ref[:, j])))
        rows.append({"trait": t, "pearson_r": r, "max_abs_diff": md})
        print("  %-5s r = %.8f   max |diff| = %.2e" % (t, r, md))
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "mtag_port_validation.csv"), index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--omega", required=True)
    ap.add_argument("--sigma", required=True)
    ap.add_argument("--mtag-dir", required=True)
    ap.add_argument("--traits", nargs="+", required=True)
    ap.add_argument("--snps", type=int, default=20000)
    ap.add_argument("--out", default="./outputs")
    a = ap.parse_args()
    main(a.omega, a.sigma, a.mtag_dir, a.traits, a.snps, a.out)
