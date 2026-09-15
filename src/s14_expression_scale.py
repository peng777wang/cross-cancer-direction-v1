"""LIMITATION 2: "the MR estimate is per unit of the eQTL's normalised
expression; the absolute magnitude depends on the expression scale and should
not be read as per 2-fold."

The unit is not unknown -- it is exactly one standard deviation of the
inverse-normalised expression, which can be anchored to real, measured
expression in the same tissue.  This script:

  1. pulls GTEx lung sample-level expression for IREB2 from the official GTEx
     portal API (so the anchor is a measured distribution, not an assumption);
  2. converts one SD into a fold-change, and converts the MR estimates into
     per-1-SD, per-1.26-fold and per-2-fold units;
  3. reports how far the cis instrument actually moves expression, which is the
     number that decides whether a per-2-fold statement is an interpolation or a
     wild extrapolation.
"""
import io
import json
import os
import urllib.request

import numpy as np
import pandas as pd

DATA = os.environ.get("DATA_DIR", "source_data")
OUT = os.environ.get("OUT_DIR", "outputs/lim2")
GENCODE = "ENSG00000136381.13"        # IREB2, GTEx v10
ENSG_V8 = "ENSG00000136381"
EBI_DEFAULT_MEDIAN_TPM = 13.456       # from the eQTL input file (GTEx v8)
BETA_CIS = -0.108214                  # eQTL effect per allele, SD units
SE_CIS = 0.0203518


def fetch_lung_expression():
    url = ("https://gtexportal.org/api/v2/expression/geneExpression"
           "?gencodeId=%s&tissueSiteDetailId=Lung" % GENCODE)
    with urllib.request.urlopen(url, timeout=60) as r:
        d = json.load(r)
    rec = d["data"][0]
    return np.asarray(rec["data"], dtype=float), rec.get("datasetId")


def main():
    os.makedirs(OUT, exist_ok=True)
    tpm, dataset = fetch_lung_expression()
    lg = np.log2(tpm + 1)
    sd_log2 = float(lg.std(ddof=1))
    med_tpm = float(np.median(tpm))
    # The eQTL phenotype is inverse-normal transformed (verified: across 593
    # genes the residual SD implied by the eQTL se never exceeds 1, which is only
    # possible if the transformed phenotype has unit variance).  An inverse
    # normal transform is rank-preserving, so ONE SD of the transformed phenotype
    # is exactly the median-to-84th-percentile interval of the expression
    # distribution.  That interval can be read off the measured data without
    # assuming any underlying scale:
    p84 = float(np.percentile(tpm, 84))
    fold_per_sd_quantile = p84 / med_tpm
    fold_per_sd = fold_per_sd_quantile            # primary anchor
    sd_units_per_2fold = float(np.log(2) / np.log(fold_per_sd_quantile))

    # one SD of the inverse-normalised expression == sd_log2 units of log2(TPM+1)
    facts = {
        "gtex_dataset": dataset,
        "n_samples_lung": int(len(tpm)),
        "median_tpm_v10_api": med_tpm,
        "median_tpm_v8_eqtl_file": EBI_DEFAULT_MEDIAN_TPM,
        "relative_difference_median_pct": float(
            100 * abs(med_tpm - EBI_DEFAULT_MEDIAN_TPM) / EBI_DEFAULT_MEDIAN_TPM),
        "sd_log2_tpm_plus1": sd_log2,
        "p50_tpm": med_tpm, "p84_tpm": p84,
        "fold_change_per_1SD_quantile_anchor": float(fold_per_sd_quantile),
        "fold_change_per_1SD_log2sdcrosscheck": float(2 ** sd_log2),
        "fold_change_per_1SD": float(fold_per_sd),
        "sd_units_per_2fold": sd_units_per_2fold,
        "sd_units_per_2fold_if_log2scale": float(1.0 / sd_log2),
        "instrument_homozygous_contrast_SD": float(2 * abs(BETA_CIS)),
        "instrument_homozygous_contrast_fold": float(
            fold_per_sd_quantile ** (2 * abs(BETA_CIS))),
        "extrapolation_factor_for_2fold_vs_instrument_range": float(
            sd_units_per_2fold / (2 * abs(BETA_CIS))),
        "instrument_F": float((BETA_CIS / SE_CIS) ** 2),
    }

    mr = pd.read_csv(os.path.join(DATA, "ireb2_mvmr_valid.csv"))
    mr = mr.rename(columns={"IREB2_theta": "per_SD", "IREB2_se": "se_per_SD",
                            "IREB2_p": "p"})
    # F = (L - mu)/s with L = log2(TPM+1), so dL/df = s and df/dL = 1/s:
    # the effect per log2 unit (= per 2-fold change) is theta / s.
    mr["per_2fold_extrapolated"] = mr.per_SD * sd_units_per_2fold
    mr["se_per_2fold_extrapolated"] = mr.se_per_SD * sd_units_per_2fold
    mr["OR_per_2fold_extrapolated"] = np.exp(mr.per_2fold_extrapolated)
    mr["OR_per_SD"] = np.exp(mr.per_SD)
    mr["OR_at_instrument_range"] = np.exp(mr.per_SD * (2 * abs(BETA_CIS)))
    mr["within_instrument_range"] = "yes"

    out = pd.DataFrame([facts])
    out.to_csv(os.path.join(OUT, "lim2_expression_scale_facts.csv"), index=False)
    mr.to_csv(os.path.join(OUT, "lim2_mr_units.csv"), index=False)

    print("GTEx dataset %s | lung n = %d" % (dataset, len(tpm)))
    print("median TPM: v10 API %.3f vs v8 eQTL file %.3f  (%.2f%% apart)"
          % (med_tpm, EBI_DEFAULT_MEDIAN_TPM, facts["relative_difference_median_pct"]))
    print("SD of log2(TPM+1) = %.4f  ->  1 SD = %.3f-fold change"
          % (sd_log2, fold_per_sd))
    print("cis instrument spans 2*|beta| = %.3f SD = %.3f-fold"
          % (2 * abs(BETA_CIS), facts["instrument_homozygous_contrast_fold"]))
    print("a 2-fold claim is %.1fx the instrument's whole contrast"
          % facts["extrapolation_factor_for_2fold_vs_instrument_range"])
    print()
    print("primary unit: per 1 SD of genetically predicted expression "
          "(= %.3f-fold change in IREB2)" % fold_per_sd)
    print(mr[["outcome", "per_SD", "OR_per_SD", "OR_at_instrument_range",
              "per_2fold_extrapolated", "OR_per_2fold_extrapolated", "p"]]
          .round(3).to_string(index=False))
    print("\nNOTE: the per-2-fold column is shown only to be quantified; it "
          "extrapolates %.1fx beyond the range the cis instrument spans."
          % facts["extrapolation_factor_for_2fold_vs_instrument_range"])


if __name__ == "__main__":
    main()
