# Audit log: every error found during the project, and how it was found

This log is part of the deliverable. Each entry states the error, how it was
detected, what it invalidated, and what replaced it. Five of the seven errors
changed a published-in-draft conclusion.

| # | error | detected by | invalidated | replaced by |
|---|---|---|---|---|
| 1 | **coloc H3 prior used `p12` instead of `p1*p2`** in a hand-written port | line-by-line comparison against `coloc::combine.abf` (`R/claudia.R` 109–133) | the claim "no colocalisation anywhere" (a false negative: H3 inflated 10,000× at default priors) | `src/s05_coloc_core.py`; re-run gave PP.H4 = 0.72–0.996 at the IREB2 locus |
| 2 | **PLINK `.bed` dosage direction assumed, not verified** | comparing dosage-derived allele frequencies with the GWAS `FRQ` column (median \|Δf\| = 0.013 vs 0.51 for the wrong convention) | potential global sign error in the LD matrix | `common.read_plink_bed_snps`, convention documented; independently re-checked against published LD values (r = 0.994, 0.957) |
| 3 | **strand not handled** in allele comparison (raw string match) | allele-frequency QC flagged 53% of SNPs as discordant | 1,555 of 2,941 SNPs were being silently dropped | complement-aware `common.harmonise`; only 0.3% now dropped |
| 4 | **`estimate_s_rss` misread** (s≈1 assumed to mean "consistent") | reading the function source: s is the *iid* component, so s≈0 means consistent | a transient, wrong conclusion that "the LD matrix does not match the data" | corrected interpretation; the LD matrix was consistent (cor(z, Rz) = 0.95) |
| 5 | **MVMR instruments were not significance-filtered** (LD-pruned only; mean F = 2.9) | computing instrument strength F before trusting the estimate | two multivariable-MR runs ("the effect halves after smoking adjustment") | instrument set rebuilt with p < 1e-5 + LD pruning (F = 28.3) plus 14 external smoking instruments; the corrected estimate is −1.35…−2.58, i.e. *unchanged* by adjustment |
| 6 | **smoking phenotype mixed up**: FinnGen `SMOKING` is a diagnosis-based endpoint (4,271 cases), not a behavioural smoking measure | only 3 of 15 established smoking loci replicate in it | any smoking-adjusted estimate using it | GSCAN smoking initiation (GWAS Catalog `GCST007474`, up to 1.23M) |
| 7 | **LD-index bug in a sensitivity script**: the metadata table was subset first, then used to index the full LD matrix | R raised "prior variance is unreasonably large … check the LD matrix" | that sensitivity script only (the main `coloc.susie` run uses the full table and was verified correct by reading the code) | fixed; main run confirmed correct |

## Errors found by *running* the code, not by reading it

Items 1–7 were found by auditing numbers and reading source. The four below were
found only by executing the released scripts end to end. They are listed
separately because reading them was not enough, which is the argument for running
the released code before submission.

| # | defect | found by | effect | fix |
|---|---|---|---|---|
| 8 | `s06_coloc_susie.R` aborted on the first gene that had no credible-set pair: `coloc.susie` returns a list with only `nsnps` (so `$summary` is `NULL`), and the guard tested the object rather than the element | executing `s06_coloc_susie.R` on the server against the real 15q25 region | the script could never finish on real data — it would have shipped broken | guard on `is.null(cc$summary)` and on `nrow(sm) == 0`; re-run reproduces PP.H4 = 0.715 / 0.962 / 0.996 |
| 9 | `s12_mtag_simulation.py` crashed when the `--out` directory did not exist | running it with a fresh output directory | the script was unrunnable from a clean checkout | `os.makedirs(out, exist_ok=True)` |
| 10 | `verify_claims.py` "verified" the cis-instrument F by comparing the constant 28.3 with 28.3 — a check that could never fail | reading the check while asking what it actually tested | none of the reported numbers, but the verification was vacuous | `source_data/instrument_strength.csv` added; F is now re-derived from each instrument's beta and se |
| 11 | the function exposed as conditional F was a residual-regression shortcut, not Sanderson's conditional F | comparing the implementation with the published definition before writing the Methods | no reported number (the manuscript reports marginal F only) | renamed to `approximate_conditional_f` with an explicit caveat; Methods states that marginal F is what is reported |

## What was verified *correct* (and how)

| item | verification |
|---|---|
| MTAG port | 20,000 SNPs × 10 traits against a real MTAG run: r = **1.00000000**, max \|Δ\| = 7×10⁻¹⁵ |
| coloc port | four synthetic scenarios (shared variant → H4; one-sided → H1; distinct → H3; none → H0) |
| LD matrix | LD(rs16969968, rs1051730) = 0.994; LD(rs55781567, rs16969968) = 0.957 — matches published EUR values |
| allele direction | three key variants checked by hand, eQTL (per ALT) vs GWAS (per A1) |
| MVMR estimator | 200 simulated replicates: bias −0.004 / +0.008 (SE 0.04 / 0.06) |
| instrument independence | measured in 1000G EUR: max pairwise \|r\| = 0.117, none > 0.2 |
| winner's curse | instrument selected in GTEx (n = 515), independent of every cancer GWAS |
| heterogeneity impact | Cochran Q reported per outcome; the propagation into the IREB2 estimate is < 0.5% because the cis variant's smoking effect is ≈ 0 |
| coloc prior sensitivity | reported as a range (0.715 / 0.962 / 0.996) rather than a single value |
| the R colocalisation script | executed on the cluster against the real 15q25 region; its PP.H4 reproduces the audited values to three decimals (see `docs/EXECUTION_LOG.md`) |
| the estimators' test suite | `tests/test_reproducibility.py`, 22 tests, covering allele harmonisation, the four coloc scenarios, Omega recovery, no sign flip under a correct model, growing borrowing with power asymmetry, MVMR unbiasedness over 500 replicates, and the 19-claim verification |

## Two claims that were *retracted* during the project

1. **"MTAG manufactures opposite-direction loci."** A simulation study (here:
   `src/s12_mtag_simulation.py`, 7 × 3 scenarios × 20 replicates × 200,000 SNPs,
   with the MTAG port validated to machine precision) shows that under a correctly
   specified model MTAG **never flips a sign**: every opposite-direction locus it
   reports is genuinely opposite in the generating model. What MTAG does is
   *amplify* a weakly powered trait's estimate by borrowing from the strong trait
   (up to ≈13% of the estimate at ρ = −0.5 with a 20-fold power difference).
   The correct statement is therefore: at loci where the weakly powered trait has
   no signal of its own, MTAG's estimate for that trait is a **model prediction,
   not an observation** — its sign may be real, but its significance comes
   entirely from the model.

2. **"The 8q24/8p21 loci show cross-cancer antagonism that is supported by the
   data."** Raw-data verification failed: the opposing arm never reaches
   \|z\| > 5.45 in its own GWAS (median \|z\| = 1.0–1.5; 68–82% of the estimate is
   borrowed), and a well-powered related cancer (colorectal, 100,204 cases) is
   *concordant* with prostate cancer at the same variant (rs6983267, both
   p < 1e-120). The antagonism claim was withdrawn; the DSI is presented as a
   descriptive, MTAG-derived annotation rather than an established biological fact.
