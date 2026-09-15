"""LIMITATION 5, sensitivity to the credible-set definition.

The region-wide test used eQTL p < 1e-5 (52 SNPs).  This recomputes the key
statistics at 1e-6 (the strongest SNPs only) and at 1e-4 (a looser set) so the
conclusion cannot be an artefact of that threshold.
"""
import os

import numpy as np
import pandas as pd

import lim5_cis_heterogeneity as L

OUT = os.path.join(L.ROOT, "lim5_out")
THRESHOLDS = [1e-6, 1e-5, 1e-4]


def main():
    meta, R, base = L.build()
    rows = []
    for trait in L.TRAITS:
        m0 = L.add_gwas(base, trait)
        ztab = pd.read_csv(os.path.join(L.ROOT, "zs_gwas_%s_15q25.tsv" % trait),
                           sep="\t", dtype={"snp": str})
        ztab.columns = ["snp", "z"]
        zmap = dict(zip(ztab.snp, ztab.z.astype(float)))
        ci = {s: k for k, s in enumerate(meta.snp)}
        for thr in THRESHOLDS:
            d = m0[m0.pvalue.astype(float) < thr].reset_index(drop=True)
            if len(d) < 3:
                rows.append({"trait": trait, "p_e_max": thr, "n_snps": len(d)})
                continue
            i = d.i.values.astype(int)
            Rs = ((R[np.ix_(i, i)] + R[np.ix_(i, i)].T) / 2)
            b_e, se_e = d.b_e.values, d.se_e.values
            b_g, se_g = d.b_g.values, d.SE.values.astype(float)
            wald = b_g / b_e
            lead = int(np.argmax(np.abs(b_e / se_e)))
            off = Rs[~np.eye(len(i), dtype=bool)]
            rng = np.random.default_rng(L.SEED)
            q, qn, qp = L.q_and_null(b_e, se_e, b_g, se_g, Rs, rng,
                                     float(wald[lead]))
            rec = {"trait": trait, "p_e_max": thr, "n_snps": len(d),
                   "mean_abs_r": float(np.abs(off).mean()),
                   "wald_median": float(np.median(wald)),
                   "wald_iqr": float(np.percentile(wald, 75) - np.percentile(wald, 25)),
                   "wald_min": float(wald.min()), "wald_max": float(wald.max()),
                   "Q": q, "Q_null_median": qn, "Q_p": qp}
            for cond in L.COND_SNPS:
                if cond in ci and cond in zmap:
                    r_jA = R[i, ci[cond]]
                    zc = (d.z_g.values.astype(float) - r_jA * zmap[cond]) / \
                        np.sqrt(np.clip(1 - r_jA ** 2, 1e-12, None))
                    rec["cond_snp"] = cond
                    rec["z_lead_raw"] = float(d.z_g.values[lead])
                    rec["z_lead_cond"] = float(zc[lead])
                    rec["max_abs_r_to_cond"] = float(np.abs(r_jA).max())
                    break
            rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "lim5_threshold_sensitivity.csv"), index=False)
    pd.set_option("display.width", 240)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
