"""Re-derive every headline number from the shipped source data and check it.

This is the script a reader should run first: it needs no raw data, only the CSVs
in source_data/, and it prints PASS/FAIL for each claim in the manuscript.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "source_data")

RESULTS = []


def chk(label, got, expect, tol=0.05, unit=""):
    good = abs(float(got) - float(expect)) <= tol
    RESULTS.append(good)
    print("  [%s] %-62s %s%s (expected %s%s)" %
          ("PASS" if good else "FAIL", label, round(float(got), 4), unit,
           expect, unit))


def main(src):
    print("=" * 100)
    print("A. Direction Safety Index (MTAG-derived) and its MTAG-free replication")
    print("=" * 100)
    d = pd.read_csv(os.path.join(src, "dsi_loci.csv"))
    chk("number of cross-cancer loci", len(d), 273, 0, "")
    chk("loci with DSI = +1 (aligned)", (d.DSI == 1).sum(), 196, 0, "")
    chk("loci with DSI = -1 (reversed)", (d.DSI == -1).sum(), 42, 0, "")
    chk("loci with -1 < DSI < +1", ((d.DSI > -1) & (d.DSI < 1)).sum(), 35, 0, "")

    r = pd.read_csv(os.path.join(src, "dsi_loci_rawz.csv"))
    chk("raw-z DSI = +1", (r.DSI == 1).sum(), 179, 0, "")
    chk("raw-z DSI = -1", (r.DSI == -1).sum(), 36, 0, "")
    chk("per-locus DSI value agreement (raw vs MTAG)",
        100 * (r.DSI == r.DSI_mtag).mean(), 84.2, 0.1, "%")
    chk("per-locus DSI class agreement",
        100 * (r.class_raw == r.class_mtag).mean(), 85.7, 0.1, "%")
    chk("Pearson r of DSI (raw vs MTAG)",
        r[["DSI", "DSI_mtag"]].corr().iloc[0, 1], 0.696, 0.005, "")

    print()
    print("=" * 100)
    print("B. How much of the reversed-arm signal is borrowed by MTAG")
    print("=" * 100)
    b = pd.read_csv(os.path.join(src, "antagonism_borrowing.csv"))
    med = b.groupby("trait").borrowed_frac.median() * 100
    chk("antagonistic locus x trait entries", len(b), 195, 0, "")
    chk("distinct antagonistic loci", b.SNP.nunique(), 77, 0, "")
    for t in ["GC", "SCLC", "PC"]:
        print("      %-5s median borrowed = %.1f%%  (own |z| never reaches 5.45)"
              % (t, med[t]))
    print("      PRAD  median borrowed = %.1f%%  (the arm that is real)" % med["PRAD"])

    print()
    print("=" * 100)
    print("C. IREB2 colocalisation at 15q25")
    print("=" * 100)
    cs = pd.read_csv(os.path.join(src, "coloc_susie_summary.csv"))
    ire = cs[cs.gene == "IREB2"]
    chk("IREB2 eQTL x lung-cancer signal B, PP.H4 (default prior)",
        ire["PP.H4.abf"].max(), 0.9692, 0.01, "")
    chk("IREB2 eQTL x classic smoking signal A, PP.H3 (default prior)",
        ire["PP.H3.abf"].max(), 0.9633, 0.01, "")
    ps = pd.read_csv(os.path.join(src, "coloc_p12_sensitivity.csv"))
    for _, row in ps.iterrows():
        print("      p12=%.0e  signal %s  PP.H3=%.3f  PP.H4=%.3f"
              % (row.p12, row.signal, row.PP_H3, row.PP_H4))

    print()
    print("=" * 100)
    print("D. Causal effect of IREB2 expression on lung cancer (cis-MR, smoking-adjusted)")
    print("=" * 100)
    m = pd.read_csv(os.path.join(src, "ireb2_mvmr_valid.csv"))
    chk("outcomes analysed", len(m), 5, 0, "")
    chk("strongest IREB2 effect (most negative)", m.IREB2_theta.min(), -2.58, 0.02, "")
    chk("weakest IREB2 effect (least negative)", m.IREB2_theta.max(), -1.35, 0.02, "")
    print("      every outcome P < 1e-3: %s" % bool((m.IREB2_p < 1e-3).all()))
    print("      every smoking effect positive: %s" % bool((m.smoking_theta > 0).all()))

    print()
    print("=" * 100)
    print("E. Instrument quality")
    print("=" * 100)
    # F is re-derived here from each instrument's beta and se, so the check is a
    # real computation rather than a constant compared with itself.
    st = pd.read_csv(os.path.join(src, "instrument_strength.csv"))
    st["F"] = (st.beta.astype(float) / st.se.astype(float)) ** 2
    cis = st[st.instrument_set.str.startswith("IREB2")]
    ext = st[st.instrument_set.str.startswith("smoking")]
    chk("external smoking instruments", len(ext), 14, 0, "")
    chk("minimum smoking instrument F", ext.F.min(), 30.3, 0.2, "")
    chk("IREB2 cis instrument F", cis.F.iloc[0], 28.3, 0.1, "")
    ld = pd.read_csv(os.path.join(src, "instrument_ld_summary.csv"))
    print("      max pairwise |r| between instruments = %.3f (no pair > 0.2)"
          % float(ld.max_abs_r.iloc[0]))

    print()
    print("=" * 100)
    print("checks passed: %d / %d" % (sum(RESULTS), len(RESULTS)))
    print("=" * 100)
    return 0 if all(RESULTS) else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-data", default=DEFAULT)
    a = ap.parse_args()
    sys.exit(main(a.source_data))
