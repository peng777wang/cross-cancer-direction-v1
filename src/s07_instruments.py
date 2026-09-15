"""Build the instrument sets and run every instrument-quality diagnostic.

Two exposures are instrumented:
  * IREB2 expression, by the lead cis-eQTL variant (the whole cis signal is a
    single independent signal, so the region supports exactly one instrument)
  * smoking initiation, by genome-wide significant GSCAN loci that lie outside
    15q25 and are clumped (>1 Mb apart, 15q25 excluded)

Reported for each set: marginal F, conditional F, and the measured pairwise LD
between instruments (so the independence assumption is checked, not assumed).
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import ld_matrix, read_plink_bed_snps
from s09_mvmr import conditional_f, marginal_f


def clump_smoking(gscan_gz, gws_p=5e-8, window=1_000_000, exclude=(15, 78_000_000, 79_500_000),
                  max_instruments=15):
    rows = []
    with __import__("gzip").open(gscan_gz, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            try:
                pv = float(p[7])
            except (IndexError, ValueError):
                continue
            if pv < gws_p:
                rows.append((p[0], int(p[1]), p[2], p[3], p[4], float(p[8]), float(p[9]), pv))
    d = pd.DataFrame(rows, columns=["CHR", "POS", "RSID", "REF", "ALT", "BETA", "SE", "P"])
    d = d.sort_values("P")
    keep = []
    for _, r in d.iterrows():
        if r.CHR == str(exclude[0]) and exclude[1] <= r.POS <= exclude[2]:
            continue
        if all(not (r.CHR == k.CHR and abs(r.POS - k.POS) < window) for k in keep):
            keep.append(r)
        if len(keep) >= max_instruments:
            break
    return pd.DataFrame(keep).reset_index(drop=True)


def instrument_ld(panel, chrom, positions, r2_threshold=0.2):
    """Return the maximum pairwise |r| among the instruments (needs the panel)."""
    bim = pd.read_csv(panel + ".bim", sep="\t", header=None,
                      names=["chr", "snp", "cm", "bp", "a1", "a2"], dtype=str)
    want = set()
    for c, p in zip(chrom, positions):
        near = bim[(bim.chr == str(c)) & (bim.bp.astype(int) == int(p))]
        want |= set(near.snp)
    idx = np.where(bim.snp.isin(want).values)[0]
    if len(idx) < 2:
        return np.nan, []
    dos = read_plink_bed_snps(panel, idx)
    R, keep = ld_matrix(dos)
    names = bim.snp.values[idx][keep]
    off = R[~np.eye(len(names), dtype=bool)]
    return float(np.abs(off).max()), list(zip(names, np.abs(off)))


def main(gscan_gz, panel, chroms, positions, max_r=0.2):
    sm = clump_smoking(gscan_gz)
    print("smoking instruments (clumped, 15q25 excluded): %d" % len(sm))
    print(sm[["RSID", "CHR", "POS", "BETA", "SE", "P"]].to_string(index=False))
    mf = marginal_f(sm.BETA.values, sm.SE.values)
    print("\nsmoking instrument strength: mean F = %.1f, min F = %.1f" % (mf, np.min(
        (sm.BETA.values / sm.SE.values) ** 2)))
    if panel:
        mx, pairs = instrument_ld(panel, chroms, positions)
        print("measured max pairwise |r| between instruments: %.3f" % mx)
        print("pairs above %.2f: %d" % (max_r, sum(1 for _, v in pairs if v > max_r)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gscan", required=True, help="GSCAN smoking-initiation .txt.gz")
    ap.add_argument("--panel", default=None, help="PLINK prefix for the LD check")
    ap.add_argument("--chrom", nargs="*", default=[])
    ap.add_argument("--pos", nargs="*", type=int, default=[])
    a = ap.parse_args()
    main(a.gscan, a.panel, a.chrom, a.pos)
