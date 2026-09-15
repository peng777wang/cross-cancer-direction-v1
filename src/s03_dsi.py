"""Direction Safety Index: locus-level sign agreement, MTAG and MTAG-free.

Input : the MTAG output files (one per trait) and the pleiotropic-locus list.
Output: dsi_loci.csv (MTAG directions) and dsi_loci_rawz.csv (raw directions),
        plus the global null tests reported in the manuscript.
"""
import argparse
import itertools
import os

import numpy as np
import pandas as pd

TRAITS = ["LUAD", "LUSC", "SCLC", "EC", "EAC", "THCA", "RCC", "PC", "GC", "PRAD"]
GWS_Z = 5.451


def dsi_from(matrix, mask):
    """matrix: (loci x traits) signed effects; mask: which traits count per locus."""
    rows = []
    for i, snp in enumerate(matrix.index):
        j = np.where(mask[i])[0]
        if len(j) < 2:
            continue
        z = matrix.values[i, j]
        prod = np.array([z[a] * z[b] for a, b in itertools.combinations(range(len(j)), 2)])
        nc, na = int((prod > 0).sum()), int((prod < 0).sum())
        rows.append({"SNP": snp, "k": len(j), "n_conc": nc, "n_ant": na,
                     "DSI": (nc - na) / (nc + na),
                     "DSI_w": prod.sum() / np.abs(prod).sum()})
    return pd.DataFrame(rows).set_index("SNP")


def poisson_binomial_tail(ps, k):
    """P(X >= k) for X = sum of independent Bernoulli(ps)."""
    d = np.zeros(len(ps) + 1)
    d[0] = 1.0
    for p in ps:
        d[1:] = d[1:] * (1 - p) + d[:-1] * p
        d[0] *= 1 - p
    return float(d[k:].sum())


def main(mtag_dir, out_dir, n_perm=5000, seed=20260913):
    z_mtag, z_raw, ns = {}, {}, {}
    for i, t in enumerate(TRAITS, start=1):
        d = pd.read_csv(os.path.join(mtag_dir, "mtag_all_trait_%d.txt" % i), sep="\t",
                        usecols=["SNP", "Z", "N", "mtag_z", "mtag_pval"], low_memory=False)
        d = d.drop_duplicates("SNP").set_index("SNP")
        z_mtag[t], z_raw[t], ns[t] = d["mtag_z"], d["Z"], d["N"]

    sig_mtag = pd.DataFrame({t: z_mtag[t].abs() > GWS_Z for t in TRAITS})
    pleio = sig_mtag.sum(axis=1)
    keep = pleio[pleio >= 2].index
    M = pd.DataFrame({t: z_mtag[t].reindex(keep) for t in TRAITS})
    R = pd.DataFrame({t: z_raw[t].reindex(keep) for t in TRAITS})
    mask = pd.DataFrame({t: M[t].abs() > GWS_Z for t in TRAITS}).values

    mtag_t = dsi_from(M, mask)
    raw_t = dsi_from(R, mask)
    out = raw_t.join(mtag_t[["DSI"]].rename(columns={"DSI": "DSI_mtag"}))
    out["class_raw"] = np.where(out.DSI == 1, "safe", np.where(out.DSI == -1, "risky", "mixed"))
    out["class_mtag"] = np.where(out.DSI_mtag == 1, "safe",
                                 np.where(out.DSI_mtag == -1, "risky", "mixed"))
    os.makedirs(out_dir, exist_ok=True)
    mtag_t.to_csv(os.path.join(out_dir, "dsi_loci.csv"))
    out.to_csv(os.path.join(out_dir, "dsi_loci_rawz.csv"))

    k = mtag_t["k"].values
    p_safe = 2.0 ** (-k)
    n_safe = int((mtag_t.DSI == 1).sum())
    n_risky = int((mtag_t.DSI == -1).sum())
    print("loci: %d | safe: %d (expected %.1f under H0, P = %.2e) | risky: %d (expected %.1f)"
          % (len(mtag_t), n_safe, p_safe.sum(), poisson_binomial_tail(p_safe, n_safe),
             n_risky, p_safe.sum()))
    print("raw-z re-computation: value agreement %.1f%% | class agreement %.1f%%"
          % (100 * (out.DSI == out.DSI_mtag).mean(),
             100 * (out.class_raw == out.class_mtag).mean()))

    # sign-flip null on the raw directions
    rng = np.random.default_rng(seed)
    Z = R.values.astype(float)
    def frac(Zm):
        c = a = 0
        for i in range(Zm.shape[0]):
            j = np.where(mask[i])[0]
            if len(j) < 2:
                continue
            z = Zm[i, j]
            prod = np.array([z[x] * z[y] for x, y in itertools.combinations(range(len(j)), 2)])
            c += int((prod > 0).sum()); a += int((prod < 0).sum())
        return c / (c + a)
    obs = frac(Z)
    null = np.array([frac(Z * rng.choice([-1.0, 1.0], size=Z.shape)) for _ in range(n_perm)])
    print("concordant-pair fraction: observed %.4f | null %.4f (sd %.4f) | Monte-Carlo P = %.4g"
          % (obs, null.mean(), null.std(), (np.sum(null >= obs) + 1) / (n_perm + 1)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtag-dir", required=True)
    ap.add_argument("--out", default="./outputs")
    ap.add_argument("--permutations", type=int, default=5000)
    a = ap.parse_args()
    main(a.mtag_dir, a.out, a.permutations)
