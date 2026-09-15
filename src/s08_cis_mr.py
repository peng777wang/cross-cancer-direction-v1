"""Single-instrument cis-MR (Wald ratio) with an instrument-strength check.

The instrument must be selected in an independent sample from the outcome GWAS;
here the cis-eQTL comes from GTEx and the outcomes are cancer GWAS, so there is no
winner's curse.
"""
import argparse

from s09_mvmr import wald_ratio


def main(beta_e, se_e, outcomes, label_e="exposure"):
    F = (beta_e / se_e) ** 2
    print("instrument: beta = %+.4f, se = %.4f -> F = %.1f" % (beta_e, se_e, F))
    if F < 10:
        print("  WARNING: F < 10, the estimate may be weak-instrument biased")
    for name, b_out, se_out in outcomes:
        th, se, p = wald_ratio(b_out, se_out, beta_e, se_e)
        print("  %-24s causal effect %+.3f (se %.3f, P = %.2e)" % (name, th, se, p))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-eqtl", type=float, required=True)
    ap.add_argument("--se-eqtl", type=float, required=True)
    ap.add_argument("--outcome", nargs=3, action="append", metavar=("NAME", "BETA", "SE"),
                    required=True)
    a = ap.parse_args()
    outs = [(n, float(b), float(s)) for n, b, s in a.outcome]
    main(a.beta_eqtl, a.se_eqtl, outs)
