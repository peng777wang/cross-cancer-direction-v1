# Direction-aware analysis of cross-cancer genetic pleiotropy

Code accompanying the manuscript. It contains everything needed to reproduce the
numbers reported in the paper **from public summary statistics**, in the order in
which the analyses were run.

## What is here

```
src/            analysis scripts (see the pipeline table below)
tests/          self-contained regression tests (stdlib unittest; no raw data)
figures/        the manuscript figure script (reads only source_data/)
resource/       the searchable cross-cancer direction + MTAG-borrowing explorer
docs/           methods notes, audit log and execution log
source_data/    one CSV per figure panel / table cell
config.example.yaml   copy to config.yaml and fill in your own paths
```

## Verify it in two commands

```bash
python tests/test_reproducibility.py     # 22 tests: estimators + all 19 claims
python src/verify_claims.py              # the 19 headline numbers, re-derived
```

Neither needs raw data. `docs/EXECUTION_LOG.md` records what each script printed
when it was actually run, on which versions, including the R script that was
executed on the cluster against the real 15q25 LD matrix.

## Pipeline

| step | script | input | output |
|---|---|---|---|
| 1 | `src/s01_build_ld.py` | 1000G EUR PLINK panel | LD matrix + allele-QC report for a region |
| 2 | `src/s02_extract_eqtl.py` | GTEx v8 cis-eQTL (tabix) | region z-scores, harmonised to the reference alleles |
| 3 | `src/s03_dsi.py` | MTAG output for 10 cancers | locus-level DSI, with the MTAG-free re-computation |
| 4 | `src/s04_borrowing.py` | MTAG output + raw z | per-locus borrowing fraction |
| 5 | `src/s05_coloc_core.py` | — | exact port of `coloc::combine.abf` (used by step 6 and by the eQTL scan) |
| 6 | `src/s06_coloc_susie.R` | LD matrix + z-scores | `coloc.susie` credible sets and PP.H4, over three p12 priors |
| 7 | `src/s07_instruments.py` | GTEx eQTL + GSCAN smoking GWAS | instrument sets (cis + external smoking loci) with strength diagnostics |
| 8 | `src/s08_cis_mr.py` | instruments + outcome GWAS | single-instrument cis-MR (Wald ratio) |
| 9 | `src/s09_mvmr.py` | instruments + outcome GWAS | multivariable MR + diagnostics (conditional F, Cochran Q) |
| 10 | `src/s10_mtag_core.py`, `src/s11_mtag_validate.py` | official MTAG output | port of the MTAG estimator, validated to machine precision |
| 11 | `src/s12_mtag_simulation.py` | — | simulation study of MTAG borrowing under power asymmetry |
| 12 | `src/s13_coloc_prior_grid.R` | LD matrix + z-scores | every coloc posterior over a 5 x 5 prior grid, plus the exact prior ratio at which H3 and H4 cross |
| 13 | `src/s14_expression_scale.py` | `source_data/ireb2_mvmr_valid.csv` + GTEx API | the eQTL SD anchored to measured expression, and the MR estimates in interpretable units |
| 14 | `src/s15_assumption_bounds.py` | `source_data/mvmr2_instruments.csv`, `own_gwas_mvmr2.tsv` | the smoking-endpoint positive control and the bias bound for unmeasured trans effects |
| 15 | `src/s16_cis_heterogeneity.py` | region LD matrix + z-scores (server-side) | sentinel independence, homogeneity across the cis credible set and the estimate conditioned on the second signal |
| 16 | `src/s17_cis_threshold_sensitivity.py` | same, with the credible-set threshold varied | the same statistics at eQTL p < 1e-6 / 1e-5 / 1e-4 |

Steps 12-16 are the analyses that resolve the five statistical limitations listed
in the manuscript; every number they produce is written up in
`docs/LIMITATIONS_RESOLVED.md`, and `docs/AUDIT_SELFCHECK.md` records the
step-by-step scientific check of those analyses, including three substantive
errors found there and corrected (an unverified effect-size scale, a p-value-only
criterion that overstated robustness, and a positive control confounded by power).

## Regenerating Fig. 1

```bash
python figures/fig1_dsi.py            # reads ../source_data, writes ./figures_out
```

It needs only `numpy`, `pandas` and `matplotlib`; it does not read the upstream
pipeline output. The script that produced the submitted figure is this one, and
running it reproduces the submitted PNG byte for byte. The render-time
panel-alignment gate is used when the helper is on `PYTHONPATH` and is otherwise
skipped with a notice.

## The searchable resource

`resource/pleiotropy_explorer.html` is a single self-contained page (340 KB, no
dependencies, works offline) that lets a reader query all 273 loci by rsID, gene,
region or cytogenetic band, filter by direction class and by the number of
cancers, and open any locus to compare its MTAG z against that cancer's own GWAS
z with the borrowed share of each estimate. `resource/README.md` documents it,
including how to host it and what it does not claim.

```bash
python ../../work/build_resource.py          # rebuild the data
python resource/build_explorer.py            # inline it into the page
node   resource/test_explorer.js             # 21 functional checks
```

## Requirements

```
python -m pip install -r requirements.txt
Rscript -e 'install.packages(c("susieR","coloc"))'
```

`config.yaml` holds every path (PLINK panel, eQTL files, GWAS files, output
directory). **No credentials and no absolute paths are hard-coded in the
scripts**, and none of the shipped files contains a password: the scripts here
never connect to a server.

## Quick verification

```bash
python src/verify_claims.py --config config.yaml
```

This re-derives every headline number from the saved source data and prints
PASS/FAIL for each, so a reader can confirm the reported values without re-running
the (very large) upstream analyses.

## Data sources (all public)

| resource | version | used for |
|---|---|---|
| GTEx v8 cis-eQTL (EBI eQTL Catalogue) | lung, prostate, stomach | eQTL instruments, colocalisation, tissue specificity |
| 1000 Genomes EUR | 503 samples | LD matrices, allele harmonisation |
| GWAS Catalog summary statistics | see `docs/METHODS_NOTES.md` | ten cancer GWAS, GSCAN smoking, FinnGen R11 |
| Open Targets Platform v4 | — | druggability, known drugs |
| DGIdb | — | approved-drug map |
| gnomAD | GRCh38 | gene constraint |

## Notes on rigour

`docs/AUDIT_LOG.md` records every error that was found and corrected during the
project (including five that invalidated earlier results), how each was detected,
and what replaced it — including four defects that were only found by *running*
the released code rather than reading it. The scripts here are the corrected
versions.

Two claims were retracted during the project and are **not** made anywhere in
this code or in the figures: that MTAG manufactures opposite-direction loci, and
that the 8q24/8p21 loci show cross-cancer antagonism supported by their own
GWAS. `docs/AUDIT_LOG.md` states what replaced each.

## Source data index

| file | panel / claim it supports |
|---|---|
| `dsi_loci.csv`, `dsi_loci_rawz.csv` | Fig. 1a, b, e — DSI classes under MTAG and under raw z |
| `antagonism_borrowing.csv` | Fig. 1c — borrowing fraction (this is what explains the reversed class) |
| `dsi_genes.csv`, `dsi_drug_targets_coding.csv` | Fig. 1d — gene-level DSI and existing drug annotations |
| `coloc_susie_summary.csv`, `coloc_p12_sensitivity.csv` | IREB2 colocalisation and its prior sensitivity |
| `coloc_all_corrected_annotated.csv` | the full eQTL x cancer colocalisation screen |
| `mvmr2_instruments.csv`, `instrument_strength.csv`, `instrument_ld_summary.csv`, `ireb2_mvmr_valid.csv` | cis-MR and smoking-adjusted MVMR with instrument diagnostics |
| `mtag_simulation_summary.csv` | the simulation behind the retraction of the "MTAG manufactures antagonism" claim |
| `eqtl_samevariant_direction_summary.csv` | cross-tissue direction concordance (99.84%) |
| `dsi_loci_positions.csv` | GRCh37 position of each DSI locus (Fig. 1a) |
| `dsi_signflip_null_1000reps.csv` | the sign-flip null behind the grey histogram in Fig. 1b (273,000 values = 273 loci x 1,000 replicates) |
| `fig1_panelF_rg.csv` | the genetic-correlation forest in Fig. 1f |
| `lim1_prior_grid.csv`, `lim1_tipping_points.csv` | the coloc prior sensitivity grid and the exact H3/H4 crossing point |
| `lim2_expression_scale_facts.csv`, `lim2_mr_units.csv` | the expression-scale anchor (1 SD = 1.226-fold, from the median-to-P84 interval) and the MR estimates in each unit |
| `lim3_trans_tipping_points.csv` | how large an unmeasured trans effect would have to be to change the conclusion |
| `lim4_positive_control_finngen_smoking.csv`, `lim4_split_half.csv` | the smoking-endpoint positive control and split-half robustness |
| `lim5_results.csv` | cis-region structure, sentinel independence, homogeneity and conditioning |
| `own_gwas_mvmr2.tsv` | the project GWAS effects for the 15 instruments |
