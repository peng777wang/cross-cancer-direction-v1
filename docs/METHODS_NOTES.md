# Methods notes: every parameter and where it comes from

This file exists so that no number in the manuscript is unexplained. Each entry
gives the value, the source or justification, and the script that uses it.

## Thresholds

| parameter | value | provenance |
|---|---|---|
| genome-wide significance | \|z\| > **5.451** | two-sided p < 5e-8: z = Φ⁻¹(1 − 2.5e-8) |
| eQTL instrument inclusion | **p < 1e-5** | region-wide Bonferroni over the 2,970 tested SNPs is 1.7e-5; 1e-5 is marginally stricter and gave a strong instrument set (F = 28.3). A variant set without this filter had mean F = 2.9 and was discarded (docs/AUDIT_LOG.md, item 5) |
| LD pruning between instruments | **r² < 0.1** | conventional independence threshold; applied greedily from the strongest eQTL variant |
| MAF filter | **> 0.01** (LD), **> 0.05** (instruments) | standard; 0.05 keeps the cis-MR instrument well imputed |
| allele-frequency concordance | **\|Δf\| < 0.15** | between 1000G EUR and the GWAS; removes multiallelic/strand-ambiguous sites |
| coloc prior p1 = p2 | **1e-4** | `coloc.abf` default |
| coloc prior p12 | **1e-5** (sensitivity over 1e-6, 1e-5, 1e-4) | `coloc.abf` default; sensitivity reported because PP.H4 depends on it |
| coloc prior sd | **0.15** quantitative, **0.2** case/control | coloc defaults (`approx.bf.estimates`) |
| susie L | **5** | number of effects allowed; the locus carries two signals |

## Data provenance

| dataset | accession / version | role |
|---|---|---|
| GTEx v8 lung cis-eQTL | EBI eQTL Catalogue, `Lung.tsv.gz` (n = 515) | exposure instrument, colocalisation |
| GTEx v8 prostate / stomach | same catalogue | tissue specificity |
| 1000 Genomes EUR | 503 samples, `EUR.bed/bim/fam` | LD reference (GRCh37) |
| LUAD / LUSC / SCLC GWAS | the project's munged summaries (see manuscript Table S1) | outcomes |
| GSCAN smoking initiation | GWAS Catalog `GCST007474` (up to 1,232,091 EUR) | second exposure (smoking) |
| FinnGen R11 | `C3_BRONCHUS_LUNG_EXALLC`, `J10_COPD`, `E4_IRON_MET`, `D3_ANAEMIA*` | surrogate outcomes, on-target safety |
| Open Targets Platform v4 | live API | tractability, known drugs |
| DGIdb | `interactions.tsv` | approved-drug map |
| gnomAD | GRCh38 | gene constraint |

## Statistical definitions

**DSI (Direction Safety Index).** For a locus significant in k cancers, take the
signs of the effects in those k cancers and form all C(k,2) pairs:

    DSI = (n_concordant - n_antagonistic) / (n_concordant + n_antagonistic) ∈ [-1, +1]

Null (independent signs): P(all k equal) = 2^(1−k); the exact Poisson-binomial
tail is used for the global test. Two sign-flip nulls are used, and they are not
the same object:

* the **per-locus DSI distribution** (grey histogram) is built from 1,000
  sign-flip replicates of the whole effect matrix, giving 273 x 1,000 = 273,000
  null DSI values;
* the **concordant-pair fraction** (448/572 under MTAG, 70.6% under raw z) is
  compared with 5,000 sign flips, so its Monte-Carlo P has a floor of 1/5,001.

Quoting a single permutation count for both would misstate one of them.

**MTAG-free re-computation.** The same DSI is recomputed with the *raw* GWAS
z-scores (the `Z` column of the MTAG input) on the same locus and trait sets, so
the comparison isolates MTAG's shrinkage.

**Borrowing fraction.** 1 − |z_raw| / |z_MTAG| for a given locus × trait. 68–82%
for the weakly powered cancers at antagonistic loci, −9% for the well-powered one.

**Colocalisation.** `susie_rss` (L = 5) per dataset, then `coloc.susie`; the
critical implementation detail is that H3 uses the prior p1·p2 and H4 uses p12
(see `src/s05_coloc_core.py`).

**cis-MR.** Wald ratio with the lead cis-eQTL variant; F = (beta/se)² = 28.27
(beta = -0.1082, se = 0.02035; re-derived in `verify_claims.py` from
`source_data/instrument_strength.csv`).
The instrument was selected in GTEx (n = 515), an independent sample from every
cancer GWAS, so selection cannot bias the outcome estimate (no winner's curse).

**Instrument strength: marginal or conditional.** The figure and the manuscript
report **marginal** F (28.3 for the cis variant; mean 55.9, minimum 30.3 across
the 14 smoking instruments). No conditional F is reported. `s09_mvmr.py` contains
an *approximate* conditional F (a residual-regression shortcut) which is named
`approximate_conditional_f` and is not the exact Sanderson (2019) statistic,
because the exact version needs the full instrument-correlation matrix; it is
included for orientation only and no reported number depends on it.

**Multivariable MR.** Two exposures (IREB2 expression, smoking initiation) with
14 external smoking instruments from GSCAN plus the cis variant; weighted least
squares with weights 1/se(Γ)². The IREB2 estimate is essentially unchanged by the
adjustment because the cis variant is not associated with smoking initiation
(p = 0.83 in GSCAN).

## Reproducing the reported numbers

| route | what it gives | cost |
|---|---|---|
| `python src/verify_claims.py` | all 19 headline numbers, re-derived from `source_data/` | seconds, no raw data |
| `python tests/test_reproducibility.py` | the above plus the estimator self-tests | ~5 s, no raw data |
| `Rscript src/s06_coloc_susie.R ...` | the colocalisation posteriors from a raw LD matrix and z-scores | minutes, needs the region inputs |
| `s01`-`s12` in order | the full pipeline from public summary statistics | hours to days, needs the raw inputs |

## Limitations, and what was done about each

These five started as caveats. Each was turned into an analysis with a stated
boundary; none is left as a bare disclaimer. Full write-up with every number:
`docs/LIMITATIONS_RESOLVED.md`.

| # | limitation | what was done | boundary now stated |
|---|---|---|---|
| 1 | PP.H4 depends on the p12 prior | evaluated a 5 x 5 prior grid; the H4/H3 odds depend on the priors only through r = p12/(p1p2), and K = (S1S2 - S12)/S12 is constant across the grid | H4 overtakes H3 at r = 38.7; PP.H4 > 0.8 needs r > 221; the coloc default r = 1000 gives PP.H4 = 0.962. For the smoking signal A, K = 5.1e4, so "distinct variants" holds throughout the conventional range |
| 2 | the MR estimate is per SD of normalised expression | verified the scale before using it: the effect is per allele (corr(se, 1/sqrt(2p(1-p))) = 0.996) and the phenotype is inverse-normal transformed (across 593 genes the implied residual SD never exceeds 1). Then anchored 1 SD using the rank-preserving median-to-P84 interval | 1 SD = 1.226-fold; the cis instrument spans only 1.045-fold; per-2-fold values are reported but extrapolate 15.7x beyond the instrumented range |
| 3 | the smoking instruments are assumed not to affect IREB2 expression | replaced the assumption with a bias bound: give all 14 instruments the same direct effect delta and find where the conclusion dies, judged by BOTH significance loss and 50% attenuation of the estimate | the direction and significance survive to delta = 0.15-0.50 SD per allele, but the MAGNITUDE halves at delta = 0.03-0.04 SD, i.e. only 28-37% of the cis effect (0.108 SD). The magnitude must therefore be reported as an upper-bound-style estimate |
| 4 | FinnGen `SMOKING` is diagnosis-based | effect-size incompatibility at the sentinel (does not depend on the endpoint's power), with the replication rate kept only as supporting evidence because the endpoint has just 4,271 cases | the same variant is beta = -0.0007 (se 0.0021) in GSCAN but +0.1400 (se 0.0258) in the endpoint: difference z = 5.43, P = 5.6e-8. No reported estimate uses the endpoint; split-half and leave-one-out leave the IREB2 estimate unchanged (< 0.01) |
| 5 | the cis region carries one independent eQTL signal | quantified the structure (52 credible-set SNPs, mean \|r\| = 0.905), the sentinel-independence of the Wald ratio, homogeneity across the credible set, a re-estimate conditioned on the independent smoking signal, and a sensitivity to the credible-set threshold (1e-6 / 1e-5 / 1e-4) | the estimate is unchanged across sentinels and thresholds (Wald median -1.80 to -1.87, -2.42, -2.55 to -2.63) and after conditioning on signal A (z = 4.07-6.04); the heterogeneity test is non-significant at every threshold (P = 0.35-0.81). The efficient multi-SNP estimator is not usable at this locus (LD effective rank 26/52), which is reported rather than claimed |

Two quantities remain unmeasured, with the reason and the boundary stated:

* the trans effects in limitation 3 (no trans-eQTL resource in the input set);
* an efficient multi-SNP estimate in limitation 5 (near-singular LD at this locus).
