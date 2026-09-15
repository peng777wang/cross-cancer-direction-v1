"""LIMITATION 5: "the cis region carries a single independent eQTL signal, so the
cis-MR relies on one variant."

Five traceable measurements, in the order they appear in the paper:

  A. STRUCTURE      how many cis SNPs carry the eQTL signal, and how correlated
                    are they?  If they are all in one LD block there is no second
                    independent instrument to obtain, which is a property of the
                    locus.
  B. SENTINEL       the Wald ratio b_out/b_exp for every credible-set SNP.  The
     INDEPENDENCE   limitation only bites if the answer depends on which variant
                    is used as the instrument.
  C. MULTI-SNP      the efficient LD-aware estimator
                        theta = (b_e' R^+ b_g)/(b_e' R^+ b_e)
                    is computed and its numerical behaviour is reported honestly:
                    with a near-singular LD matrix it is unstable, and that
                    instability is a property of the locus, not a result.
  D. RESIDUAL       the region contains a second, independent GWAS signal (the
     HOMOGENEITY    smoking signal A).  Heterogeneity across cis SNPs is
                    therefore expected; the test is repeated on GWAS effects
                    *conditioned on signal A* (GCTA-style single-SNP conditioning).
                    If homogeneity is restored, the residual cis-MR is a clean
                    single-causal-variant estimate.
  E. LEAVE-ONE-OUT  drop the sentinel, drop the top five.

`b(e|g)` are harmonised to the LD panel's A1, so all LD correlations and all
effect sizes share one orientation.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.environ.get("CGS_ROOT", "data")   # path to the region inputs
OUT = os.path.join(ROOT, "lim5_out")
os.makedirs(OUT, exist_ok=True)

TRAITS = ["LUAD", "LUSC", "SCLC"]
WINDOW = (78_500_000, 79_250_000)
P_E_MAX = 1e-5
COND_SNPS = ["rs16969968", "rs4887067", "rs1051730"]   # signal A candidates
N_BOOT = 1000
N_MC = 2000
SEED = 20260914
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def same(a, b):
    a, b = str(a).upper(), str(b).upper()
    return a == b or a == COMP.get(b, "?")


def orient(effect, other, ref_a1, ref_a2):
    if same(effect, ref_a1) and same(other, ref_a2):
        return 1
    if same(effect, ref_a2) and same(other, ref_a1):
        return -1
    return 0


def build():
    meta = pd.read_csv(os.path.join(ROOT, "ld_15q25_meta.tsv"), sep="\t",
                       dtype={"snp": str})
    R = pd.read_csv(os.path.join(ROOT, "ld_15q25_R.tsv"), sep="\t", index_col=0)
    assert list(R.index) == list(meta.snp), "LD matrix and metadata are misaligned"
    R = R.values.astype(float)
    eq = pd.read_csv(os.path.join(ROOT, "ireb2_lung_eqtl_region.tsv"), sep="\t",
                     dtype={"rsid": str})
    eq = eq[eq.gene_id == "ENSG00000136381"].dropna(subset=["beta", "se", "rsid"])
    eq = eq.drop_duplicates("rsid")
    idx = {s: i for i, s in enumerate(meta.snp)}
    m = meta[["snp", "a1", "a2", "bp"]].merge(
        eq[["rsid", "beta", "se", "ref", "alt", "pvalue"]], left_on="snp",
        right_on="rsid", how="inner")
    m = m[(m.bp >= WINDOW[0]) & (m.bp <= WINDOW[1]) & m.snp.isin(idx)].copy()
    m = m.assign(o=[orient(r.alt, r.ref, r.a1, r.a2) for r in m.itertuples()])
    m = m[m.o != 0].copy()
    m["b_e"] = m.beta.astype(float) * m.o
    m["se_e"] = m.se.astype(float)
    m["i"] = m.snp.map(idx)
    return meta, R, m


def add_gwas(m, trait):
    z = pd.read_csv(os.path.join(ROOT, "zs_gwas_%s_15q25.tsv" % trait), sep="\t",
                    dtype={"snp": str})
    z.columns = ["snp", "z_g"]
    m = m.merge(z, on="snp", how="left")
    g = pd.read_csv(os.path.join(ROOT, "clean", "%s.txt.gz" % trait), sep="\t",
                    usecols=["SNP", "CHR", "POS", "SE"], dtype={"SNP": str})
    g = g[(g.CHR == 15) & (g.POS >= WINDOW[0]) & (g.POS <= WINDOW[1])]
    g = g.drop_duplicates("SNP")[["SNP", "SE"]].rename(columns={"SNP": "snp"})
    m = m.merge(g, on="snp", how="left")
    m = m.dropna(subset=["z_g", "SE"])
    m["b_g"] = m.z_g.astype(float) * m.SE.astype(float)
    return m


def theta_smr(b_e, b_g, R, rcond=1e-6):
    u = np.linalg.pinv(R, rcond=rcond) @ b_e
    return float(u @ b_g / (u @ b_e))


def q_and_null(b_e, se_e, b_g, se_g, R, rng, theta_fixed, n_mc=N_MC):
    """Cochran Q across correlated cis SNPs + LD-aware Monte-Carlo null.

    theta is held FIXED at a prespecified value (the sentinel's Wald ratio) rather
    than re-estimated inside the statistic.  Re-estimating it would use the
    pseudo-inverse of a numerically singular LD matrix, which makes the null
    distribution itself diverge (documented in the output as `theta_smr_boot_sd`).
    """
    w = 1.0 / (se_g ** 2 + theta_fixed ** 2 * se_e ** 2)
    q_obs = float(np.sum(w * (b_g - theta_fixed * b_e) ** 2))
    cov = np.outer(se_g, se_g) * R
    ev, V = np.linalg.eigh((cov + cov.T) / 2)
    L = V * np.sqrt(np.clip(ev, 0, None))
    qs = np.empty(n_mc)
    for k in range(n_mc):
        bg = theta_fixed * b_e + L @ rng.standard_normal(len(b_g))
        qs[k] = float(np.sum(w * (bg - theta_fixed * b_e) ** 2))
    return q_obs, float(np.median(qs)), float((np.sum(qs >= q_obs) + 1) / (n_mc + 1))


def run_trait(base, meta, R, trait, rng):
    m = add_gwas(base, trait)
    d = m[m.pvalue.astype(float) < P_E_MAX].reset_index(drop=True)
    i = d.i.values.astype(int)
    Rs = ((R[np.ix_(i, i)] + R[np.ix_(i, i)].T) / 2)
    b_e, se_e = d.b_e.values, d.se_e.values
    b_g, se_g = d.b_g.values, d.SE.values.astype(float)
    z_g = d.z_g.values.astype(float)
    off = Rs[~np.eye(len(i), dtype=bool)]

    res = {"trait": trait, "n_cis_tested": int(len(m)), "n_credible_set": int(len(d)),
           "mean_abs_r": float(np.abs(off).mean()), "min_abs_r": float(np.abs(off).min()),
           "max_abs_r": float(np.abs(off).max())}
    for thr in (1e-6, 1e-5, 1e-4, 1e-3):
        res["n_p_e_lt_%.0e" % thr] = int((m.pvalue.astype(float) < thr).sum())

    # ---- B. sentinel independence
    wald = b_g / b_e
    wt = pd.DataFrame({"snp": d.snp.values, "bp": d.bp.values,
                       "p_e": d.pvalue.astype(float).values,
                       "z_e": b_e / se_e, "F": (b_e / se_e) ** 2,
                       "z_g": z_g, "wald": wald})
    wt.to_csv(os.path.join(OUT, "lim5_wald_by_snp_%s.tsv" % trait), sep="\t",
              index=False)
    lead = int(np.argmax(np.abs(b_e / se_e)))
    res.update(lead_snp=d.snp.values[lead], lead_wald=float(wald[lead]),
               lead_F=float((b_e[lead] / se_e[lead]) ** 2),
               wald_median=float(np.median(wald)),
               wald_iqr=float(np.percentile(wald, 75) - np.percentile(wald, 25)),
               wald_min=float(wald.min()), wald_max=float(wald.max()))

    # ---- C. efficient estimator: compute it, and report why it is unusable here
    ev = np.linalg.eigvalsh(Rs)
    res["ld_effective_rank"] = int((ev > 1e-6 * ev.max()).sum())
    res["ld_condition_number"] = float(ev.max() / max(ev.min(), 1e-300))
    res["theta_smr_point"] = theta_smr(b_e, b_g, Rs)
    bs = np.array([theta_smr(b_e + se_e * rng.standard_normal(len(b_e)),
                            b_g + se_g * rng.standard_normal(len(b_g)), Rs)
                   for _ in range(N_BOOT)])
    res["theta_smr_boot_sd"] = float(np.nanstd(bs, ddof=1))

    # ---- D. does conditioning on the independent smoking signal remove the
    #         heterogeneity?
    ztab = pd.read_csv(os.path.join(ROOT, "zs_gwas_%s_15q25.tsv" % trait),
                       sep="\t", dtype={"snp": str})
    ztab.columns = ["snp", "z"]
    zmap = dict(zip(ztab.snp, ztab.z.astype(float)))
    snp_index = {s: k for k, s in enumerate(meta.snp)}
    for cond in COND_SNPS:
        if cond not in snp_index or cond not in zmap:
            continue
        ci = snp_index[cond]
        r_jA = R[i, ci]
        zA = zmap[cond]
        denom = np.sqrt(np.clip(1 - r_jA ** 2, 1e-12, None))
        z_cond = (z_g - r_jA * zA) / denom
        b_gc = z_cond * se_g
        th_fix = float(wald[lead])
        q_raw, qn_raw, p_raw = q_and_null(b_e, se_e, b_g, se_g, Rs, rng, th_fix)
        q_cnd, qn_cnd, p_cnd = q_and_null(b_e, se_e, b_gc, se_g, Rs, rng, th_fix)
        res["cond_snp"] = cond
        res["Q_raw"], res["Q_null_raw"], res["Q_p_raw"] = q_raw, qn_raw, p_raw
        res["Q_cond"], res["Q_null_cond"], res["Q_p_cond"] = q_cnd, qn_cnd, p_cnd
        res["z_lead_cond"] = float(z_cond[lead])
        res["z_lead_raw"] = float(z_g[lead])
        res["rho_wald_vs_rA"] = float(np.corrcoef(wald, r_jA)[0, 1])
        res["theta_cond_median"] = float(np.median(b_gc / b_e))
        break

    # ---- E. leave-one-out
    keep = np.ones(len(d), dtype=bool); keep[lead] = False
    res["wald_median_no_sentinel"] = float(np.median(wald[keep]))
    top5 = np.argsort(-np.abs(b_e / se_e))[:5]
    keep5 = np.ones(len(d), dtype=bool); keep5[top5] = False
    res["wald_median_no_top5"] = float(np.median(wald[keep5]))
    return res, wt


def main():
    meta, R, base = build()
    rows, tabs = [], []
    for trait in TRAITS:
        rng = np.random.default_rng(SEED)
        r, wt = run_trait(base, meta, R, trait, rng)
        wt["trait"] = trait
        rows.append(r); tabs.append(wt)
        print("=" * 74)
        for k, v in r.items():
            print("  %-26s %s" % (k, v))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "lim5_results.csv"), index=False)
    pd.concat(tabs).to_csv(os.path.join(OUT, "lim5_wald_all.tsv"), sep="\t", index=False)
    print("\nwrote", os.path.join(OUT, "lim5_results.csv"))


if __name__ == "__main__":
    main()
