"""Build a 1000G-EUR LD matrix for a region and harmonise GWAS z-scores onto it.

Verification of the two things that most easily go wrong:
  * the PLINK dosage direction (see common.read_plink_bed_snps)
  * the strand of every allele (see common.harmonise)
The script prints, for a few well-known variant pairs, the measured LD so that it
can be compared with published values before the matrix is used downstream.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import harmonise, ld_matrix, read_plink_bed_snps


def main(panel, chrom, start, end, gwas_dir, traits, out_dir,
         check_pairs=("rs16969968:rs1051730",)):
    bim = pd.read_csv(panel + ".bim", sep="\t", header=None,
                      names=["chr", "snp", "cm", "bp", "a1", "a2"], dtype=str)
    sel = (bim.chr == str(chrom)) & (bim.bp.astype(int) >= start) & (bim.bp.astype(int) <= end)
    idx = np.where(sel.values)[0]
    dos = read_plink_bed_snps(panel, idx)
    R, keep = ld_matrix(dos)
    meta = bim.iloc[idx[keep]].reset_index(drop=True)
    meta["f_ref_A1"] = np.nanmean(dos[keep], axis=1) / 2
    print("reference variants in the window: %d (kept %d)" % (len(idx), len(meta)))

    loc = {s: i for i, s in enumerate(meta.snp)}
    for pair in check_pairs:
        a, b = pair.split(":")
        if a in loc and b in loc:
            print("  LD(%s, %s) = %.3f" % (a, b, R[loc[a], loc[b]]))

    os.makedirs(out_dir, exist_ok=True)
    np.save(os.path.join(out_dir, "ld_region_R.npy"), R)
    meta.to_csv(os.path.join(out_dir, "ld_region_meta.tsv"), sep="\t", index=False)

    for tr in traits:
        g = pd.read_csv(os.path.join(gwas_dir, "%s.txt.gz" % tr), sep="\t",
                        usecols=["SNP", "CHR", "POS", "A1", "A2", "BETA", "SE", "FRQ", "N"],
                        dtype={"SNP": str, "A1": str, "A2": str})
        g = g[(g.CHR == int(chrom)) & (g.POS >= start) & (g.POS <= end)].drop_duplicates("SNP")
        m = meta.merge(g, left_on="snp", right_on="SNP", how="inner")
        z, ok = [], []
        for _, r in m.iterrows():
            flip, ambiguous, good = harmonise(r.A1, r.A2, r.a1, r.a2)
            if not good:
                ok.append(False); z.append(np.nan); continue
            if ambiguous:                      # resolve with allele frequency
                d_no = abs(r.f_ref_A1 - r.FRQ)
                d_fl = abs(r.f_ref_A1 - (1 - r.FRQ))
                flip = d_fl < d_no
            if abs(r.f_ref_A1 - ((1 - r.FRQ) if flip else r.FRQ)) > 0.15:
                ok.append(False); z.append(np.nan); continue
            ok.append(True)
            z.append((r.BETA / r.SE) * (-1 if flip else 1))
        s = pd.Series(np.nan, index=meta.snp)
        s.loc[m.snp.values] = np.array(z)
        s.to_csv(os.path.join(out_dir, "z_%s_region.tsv" % tr), sep="\t",
                 index_label="snp", na_rep="NA")
        print("  %-5s harmonised %d / %d variants (n = %.0f)"
              % (tr, int(np.sum(ok)), len(m), float(g.N.median())))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, help="PLINK prefix, e.g. /path/EUR")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--gwas-dir", required=True)
    ap.add_argument("--traits", nargs="+", required=True)
    ap.add_argument("--out", default="./outputs")
    a = ap.parse_args()
    main(a.panel, a.chrom, a.start, a.end, a.gwas_dir, a.traits, a.out)
