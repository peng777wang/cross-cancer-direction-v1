"""LIMITATION 3 and LIMITATION 4.

LIMITATION 3: the 14 external smoking instruments are assumed to have no effect
on IREB2 expression, and trans effects cannot be tested with cis-only eQTL data.
The assumption cannot be measured here, so it is bounded instead: give every
smoking instrument the same direct effect delta on IREB2 expression (the worst
case for a common trans effect), recompute the multivariable IREB2 estimate over
a grid of delta, and report the delta at which the conclusion disappears.

LIMITATION 4: FinnGen's SMOKING endpoint is diagnosis-based and disagrees with
GSCAN.
  (a) POSITIVE CONTROL - the 14 instruments are genome-wide significant loci from
      a 1.2M-sample behavioural smoking GWAS; if FinnGen SMOKING measured
      smoking they should replicate in it.
  (b) NO RESULT DEPENDS ON IT - every estimate uses GSCAN. Robustness of the
      smoking term is checked by splitting the instruments in half and by
      leave-one-instrument-out.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

DATA = os.environ.get("DATA_DIR", "source_data")
OUT = os.environ.get("OUT_DIR", "outputs/lim34")
INST = os.path.join(DATA, "mvmr2_instruments.csv")
OWN = os.path.join(DATA, "own_gwas_mvmr2.tsv")
CIS = "rs7359276"
B_E_CIS, SE_E_CIS = -0.108214, 0.0203518

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def sa(a, b):
    return a == b or a == COMP.get(b, "?")


def orient(a1, a2, ref, alt):
    if sa(a1, alt) and sa(a2, ref):
        return 1
    if sa(a1, ref) and sa(a2, alt):
        return -1
    return 0


def mvmr(B, gam, sgam):
    W = np.diag(1.0 / np.asarray(sgam, float) ** 2)
    cov = np.linalg.inv(B.T @ W @ B)
    return cov @ (B.T @ W @ np.asarray(gam, float)), cov


def build():
    inst = pd.read_csv(INST).rename(columns={"BETA": "b_smk_gscan",
                                             "SE": "se_smk_gscan"})
    own = pd.read_csv(OWN, sep="\t", dtype={"SNP": str})
    return inst, own


def frame(inst, trait_effects, name):
    m = inst.merge(trait_effects, left_on="RSID", right_on="SNP")
    o = [orient(str(a).upper(), str(b).upper(), str(c).upper(), str(d).upper())
         for a, b, c, d in zip(m.A1, m.A2, m.REF, m.ALT)]
    m = m.assign(o=o)
    m = m[m.o != 0].copy()
    m["gam"] = m.BETA.astype(float) * m.o
    m["sgam"] = m.SE.astype(float)
    m["b1"] = np.where(m.RSID == CIS, B_E_CIS, 0.0)
    m["b2"] = m.b_smk_gscan.astype(float)
    m["outcome"] = name
    return m


def estimate(m, b1=None):
    B = np.column_stack([m.b1.values if b1 is None else b1, m.b2.values])
    th, cov = mvmr(B, m.gam.values, m.sgam.values)
    se = np.sqrt(np.diag(cov))
    return th, se, th / se


def main():
    os.makedirs(OUT, exist_ok=True)
    inst, own = build()
    ext = inst[inst.RSID != CIS]

    p = ext.p_smk.astype(float).values
    n05, n08 = int((p < 0.05).sum()), int((p < 5e-8).sum())
    pc = pd.DataFrame({"rsid": ext.RSID.values, "gscan_beta": ext.b_smk_gscan.values,
                       "gscan_p": ext.P.values.astype(float),
                       "finngen_smoking_beta": ext.b_smk.values,
                       "finngen_smoking_p": p, "replicates_p05": p < 0.05,
                       "replicates_p5e8": p < 5e-8})
    pc.to_csv(os.path.join(OUT, "lim4_positive_control_finngen_smoking.csv"), index=False)
    binom_p = float(stats.binom.sf(n05 - 1, len(ext), 0.8))
    print("== LIMITATION 4a: do GSCAN smoking loci replicate in FinnGen SMOKING?")
    print("   (supporting evidence only: this endpoint has 4,271 cases, so a low")
    print("    replication rate may reflect power rather than mislabelling)")
    print("   instruments tested    : %d" % len(ext))
    print("   replicate at p < 0.05 : %d   (binomial P vs 80%% expected = %.3g)"
          % (n05, binom_p))
    print("   replicate at p < 5e-8 : %d" % n08)
    print("   cis variant in that endpoint: p = %.2g (GSCAN p = %.2g)"
          % (float(inst.loc[inst.RSID == CIS, "p_smk"].iloc[0]),
             float(inst.loc[inst.RSID == CIS, "P"].iloc[0])))

    # Effect-size incompatibility at the sentinel: this argument does not depend
    # on the endpoint's power, only on the two effect estimates and their SEs.
    row = inst[inst.RSID == CIS].iloc[0]
    b_g, se_g = float(row.b_smk_gscan), float(row.se_smk_gscan)
    b_f, se_f = float(row.b_smk), float(row.se_smk)
    diff, se_diff = b_f - b_g, np.hypot(se_f, se_g)
    print("   effect-size test at %s:" % CIS)
    print("     GSCAN smoking initiation      beta %+.4f (se %.4f) -> OR %.3f "
          "[%.3f, %.3f]" % (b_g, se_g, np.exp(b_g), np.exp(b_g - 1.96 * se_g),
                            np.exp(b_g + 1.96 * se_g)))
    print("     FinnGen SMOKING endpoint      beta %+.4f (se %.4f) -> OR %.3f "
          "[%.3f, %.3f]" % (b_f, se_f, np.exp(b_f), np.exp(b_f - 1.96 * se_f),
                            np.exp(b_f + 1.96 * se_f)))
    print("     difference = %+.4f (se %.4f) -> z = %.2f, P = %.2g"
          % (diff, se_diff, diff / se_diff, 2 * stats.norm.sf(abs(diff / se_diff))))
    print("     i.e. the endpoint's effect at this variant is excluded by the")
    print("     behavioural GWAS at z = %.1f." % (diff / se_diff))

    frames = []
    for trait in ["LUAD", "LUSC", "SCLC"]:
        sub = own[own.trait == trait][["SNP", "A1", "A2", "BETA", "SE"]]
        frames.append(frame(inst, sub, "own:" + trait))
    for tag, lab in [("lungca", "FinnGen lung cancer"), ("copd", "FinnGen COPD")]:
        sub = pd.DataFrame({"SNP": inst.RSID, "A1": inst["alt_" + tag],
                            "A2": inst["ref_" + tag],
                            "BETA": inst["b_" + tag].astype(float),
                            "SE": inst["se_" + tag].astype(float)})
        frames.append(frame(inst, sub, "fg:" + lab))

    print("\n== LIMITATION 4b: is the smoking term driven by particular instruments?")
    ext_ids = list(inst.RSID.values[inst.RSID.values != CIS])
    halves = {"A": ext_ids[0::2], "B": ext_ids[1::2]}
    rows = []
    for m in frames:
        th, se, z = estimate(m)
        rec = {"outcome": m.outcome.iloc[0], "n_iv": len(m),
               "IREB2_theta": th[0],
               "IREB2_p": 2 * stats.norm.sf(abs(z[0])),
               "smoking_theta": th[1],
               "smoking_p": 2 * stats.norm.sf(abs(z[1]))}
        for h, ids in halves.items():
            mm = m[m.RSID.isin(ids + [CIS])]
            th2, _, z2 = estimate(mm)
            rec["IREB2_theta_half" + h] = th2[0]
            rec["IREB2_p_half" + h] = 2 * stats.norm.sf(abs(z2[0]))
            rec["smoking_theta_half" + h] = th2[1]
            rec["smoking_p_half" + h] = 2 * stats.norm.sf(abs(z2[1]))
        loo = []
        for r in m.RSID.values:
            if r == CIS:
                continue
            th3, _, _ = estimate(m[m.RSID != r])
            loo.append(th3[0])
        rec["IREB2_theta_loo_min"] = float(np.min(loo))
        rec["IREB2_theta_loo_max"] = float(np.max(loo))
        rows.append(rec)
    split = pd.DataFrame(rows)
    split.to_csv(os.path.join(OUT, "lim4_split_half.csv"), index=False)
    print(split.round(3).to_string(index=False))

    print("\n== LIMITATION 3: how large a direct smoking->IREB2 effect would be needed")
    deltas = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05,
              0.075, 0.1, 0.15, 0.2, 0.5, 1.0]
    br = []
    for m in frames:
        b1 = np.where(m.RSID.values == CIS, B_E_CIS, 0.0)
        base = np.where(m.RSID.values == CIS, 0.0, 1.0)
        for delta in deltas:
            th, _, z = estimate(m, b1=b1 + base * delta)
            br.append({"outcome": m.outcome.iloc[0], "delta_direct_effect_SD": delta,
                       "IREB2_theta": th[0],
                       "IREB2_p": 2 * stats.norm.sf(abs(z[0]))})
    bnd = pd.DataFrame(br)
    bnd.to_csv(os.path.join(OUT, "lim3_trans_bias_bound.csv"), index=False)
    tips = []
    for lab, g in bnd.groupby("outcome"):
        g = g.sort_values("delta_direct_effect_SD")
        th0 = float(g.IREB2_theta.iloc[0])
        # criterion 1: significance lost.  criterion 2 (the honest one for a claim
        # about magnitude): the estimate falls below half its original size.  A
        # large delta can drive theta towards 0 while the SE shrinks too, leaving a
        # spuriously "significant" near-null estimate -- hence both criteria.
        sig_lost = g[(g.IREB2_theta >= 0) | (g.IREB2_p > 0.05)]
        half_lost = g[g.IREB2_theta.abs() < 0.5 * abs(th0)]
        tips.append({"outcome": lab,
                     "delta_for_significance_loss": float(sig_lost.delta_direct_effect_SD.iloc[0])
                     if len(sig_lost) else np.inf,
                     "delta_for_half_attenuation": float(half_lost.delta_direct_effect_SD.iloc[0])
                     if len(half_lost) else np.inf,
                     "theta_at_delta_0": th0,
                     "theta_at_delta_0.05": float(g.loc[np.isclose(g.delta_direct_effect_SD, 0.05), "IREB2_theta"].iloc[0]),
                     "theta_at_delta_0.2": float(g.loc[np.isclose(g.delta_direct_effect_SD, 0.2), "IREB2_theta"].iloc[0])})
    tip = pd.DataFrame(tips)
    tip.to_csv(os.path.join(OUT, "lim3_trans_tipping_points.csv"), index=False)
    print(tip.round(3).to_string(index=False))
    print("\nreference: cis eQTL effect per allele = %.3f SD; the 14 smoking instruments"
          % abs(B_E_CIS))
    print("have GSCAN effects of %.4f to %.4f (per allele)."
          % (ext.b_smk_gscan.min(), ext.b_smk_gscan.max()))


if __name__ == "__main__":
    main()
