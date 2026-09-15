"""How much of each trait's MTAG estimate is borrowed rather than observed.

For every locus x trait entry at an antagonistic locus we compare the trait's own
GWAS z (the MTAG input `Z`) with the MTAG z. The raw z *is* that trait's own
evidence; everything beyond it is borrowed through the genetic covariance matrix.

Output: antagonism_borrowing.csv (one row per locus x trait) and a per-trait summary.
"""
import argparse
import itertools
import os

import numpy as np
import pandas as pd

TRAITS = ["LUAD", "LUSC", "SCLC", "EC", "EAC", "THCA", "RCC", "PC", "GC", "PRAD"]
GWS_Z = 5.451


def main(mtag_dir, out_dir):
    z_mtag, z_raw = {}, {}
    for i, t in enumerate(TRAITS, start=1):
        d = pd.read_csv(os.path.join(mtag_dir, "mtag_all_trait_%d.txt" % i), sep="\t",
                        usecols=["SNP", "Z", "mtag_z"], low_memory=False)
        d = d.drop_duplicates("SNP").set_index("SNP")
        z_mtag[t], z_raw[t] = d["mtag_z"], d["Z"]
    sig = pd.DataFrame({t: z_mtag[t].abs() > GWS_Z for t in TRAITS})
    keep = sig.index[sig.sum(axis=1) >= 2]
    M = pd.DataFrame({t: z_mtag[t].reindex(keep) for t in TRAITS})
    R = pd.DataFrame({t: z_raw[t].reindex(keep) for t in TRAITS})

    rows = []
    for snp in M.index:
        z = M.loc[snp].values.astype(float)
        j = np.where(np.abs(z) > GWS_Z)[0]
        if len(j) < 2 or not (z[j].max() > 0 and z[j].min() < 0):
            continue
        for i in j:
            rows.append({"SNP": snp, "trait": TRAITS[i], "z_mtag": z[i],
                         "z_raw": float(R.loc[snp, TRAITS[i]]),
                         "side": "up" if z[i] > 0 else "down"})
    d = pd.DataFrame(rows)
    d["own_significant"] = np.abs(d.z_raw) > GWS_Z
    d["sign_match"] = np.sign(d.z_raw) == np.sign(d.z_mtag)
    d["borrowed_frac"] = (np.abs(d.z_mtag) - np.abs(d.z_raw)) / np.abs(d.z_mtag)
    os.makedirs(out_dir, exist_ok=True)
    d.to_csv(os.path.join(out_dir, "antagonism_borrowing.csv"), index=False)

    g = d.groupby("trait").agg(entries=("SNP", "size"),
                               own_significant=("own_significant", "sum"),
                               median_own_z=("z_raw", lambda s: np.abs(s).median()),
                               median_mtag_z=("z_mtag", lambda s: np.abs(s).median()),
                               median_borrowed=("borrowed_frac", "median"),
                               pct_sign_match=("sign_match", lambda s: 100 * s.mean()))
    print("antagonistic entries: %d across %d loci" % (len(d), d.SNP.nunique()))
    print(g.round(3).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtag-dir", required=True)
    ap.add_argument("--out", default="./outputs")
    a = ap.parse_args()
    main(a.mtag_dir, a.out)
