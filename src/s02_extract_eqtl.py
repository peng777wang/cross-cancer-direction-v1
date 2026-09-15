"""Extract a cis-eQTL region and harmonise it onto the LD reference alleles.

Ambiguous (A/T, C/G) variants are dropped for eQTL data because the eQTL files do
not carry a per-variant allele frequency for the effect allele; keeping them
unresolved would risk silent sign errors.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import harmonise, read_eqtl_region


def main(eqtl_path, chrom, start, end, genes, meta_tsv, out_tsv):
    meta = pd.read_csv(meta_tsv, sep="\t", dtype={"snp": str})
    meta["a1"] = meta.a1.str.upper()
    meta["a2"] = meta.a2.str.upper()
    e = read_eqtl_region(eqtl_path, chrom, start, end)
    e = e[e.type == "SNP"].dropna(subset=["beta", "se", "pvalue"])
    out = pd.DataFrame(index=meta.snp)
    out.index.name = "snp"
    for ensg, sym in genes.items():
        sub = e[e.gene_id == ensg].drop_duplicates("rsid")
        if len(sub) < 100:
            print("  %-8s skipped (only %d variants)" % (sym, len(sub)))
            continue
        m = sub.merge(meta[["snp", "a1", "a2"]], left_on="rsid", right_on="snp")
        z = np.full(len(m), np.nan)
        keep = np.zeros(len(m), dtype=bool)
        for i, r in enumerate(m.itertuples()):
            flip, ambiguous, good = harmonise(r.alt, r.ref, r.a1, r.a2)
            if not good or ambiguous:
                continue
            keep[i] = True
            z[i] = (r.beta / r.se) * (-1 if flip else 1)
        s = pd.Series(np.nan, index=meta.snp)
        s.loc[m.snp.values[keep]] = z[keep]
        out[sym] = s.values
        print("  %-8s %d variants harmonised, top |z| = %.2f"
              % (sym, int(keep.sum()), np.nanmax(np.abs(z))))
    out.to_csv(out_tsv, sep="\t", na_rep="NA")
    print("wrote", out_tsv)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eqtl", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--gene", nargs=2, action="append", metavar=("ENSG", "SYMBOL"),
                    required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    main(a.eqtl, a.chrom, a.start, a.end, dict(a.gene), a.meta, a.out)
