# Execution log: what was actually run, on which versions, with which result

Reading code cannot show that it runs. Every script in `src/` was therefore
executed, and the results below are the recorded output rather than an
expectation. Anything that could not be executed in this environment is marked.

## Environment

| component | version |
|---|---|
| Python | 3.12.14 |
| numpy / pandas / scipy | 2.3.5 / 3.0.1 / 1.18.1 |
| R | 4.4.1 (2024-06-14) |
| coloc | 5.2.3 |
| susieR | 0.14.2 |
| data.table | 1.15.4 |
| LD reference | 1000G EUR, 503 samples, GRCh37 |

Python scripts need only `requirements.txt`. The R script needs `coloc` and
`susieR`; on the cluster these are read from a user library, so the script's
`Usage` line assumes they are on `.libPaths()`.

## Executed locally (no cluster, no raw data)

| what | command | result |
|---|---|---|
| coloc port self-test | `python src/s05_coloc_core.py` | shared variant -> H4 = 0.890; one-sided -> H1 maximal; two distinct variants -> H3 > H4; neither -> H0 = 0.848; sign flip leaves the posterior unchanged |
| MTAG port, method of moments | `python src/s12_mtag_simulation.py --reps 2 --snps 20000` | Omega recovered within 15% of the generating value at 20,000 SNPs |
| claim verification | `python src/verify_claims.py` | **19 / 19 PASS**, exit 0 |
| full test suite | `python tests/test_reproducibility.py` | **22 tests, OK**, covering harmonisation, the four coloc scenarios, H3-prior algebra, Omega recovery, no sign flip under a correct model, borrowing that grows with power asymmetry, MVMR unbiasedness (500 replicates, \|bias\| < 0.02 and < 3 Monte-Carlo SE), and the claim verification |

## Executed on the cluster (real data)

The R colocation script was run against the real 15q25 region with the released
inputs — a 2,941 x 2,941 LD matrix, harmonised LUAD z-scores and GTEx lung
eQTL z-scores (n = 515) — through the same command line documented in the script
header:

```
Rscript src/s06_coloc_susie.R <ld.rds> <meta.tsv> <gwas.tsv> <eqtl.tsv> <n.txt> <out.csv>
```

Output (538 SNPs in the window; 378 with both traits present; 7 genes tested):

| gene | credible pair | p12 = 1e-6 | 1e-5 | 1e-4 |
|---|---|---|---|---|
| IREB2 x signal B (`rs28498264`) | eQTL CS1 - GWAS CS2 | PP.H4 = 0.715 | **0.962** | 0.996 |
| IREB2 x signal A (`rs4887067`) | eQTL CS1 - GWAS CS1 | PP.H3 = 0.980 | 0.963 | 0.824 |
| CHRNA5 (both pairs) | - | PP.H3 > 0.98 throughout | | |
| CHRNA3, PSMA4, HYKK, MORF4L1, CTSH | no credible set in either trait | not tested | | |

This reproduces the audited numbers (0.715 / 0.962 / 0.996) and confirms that
signal A is dominated by H3, i.e. two distinct causal variants rather than one
shared variant.

The first run of this script **failed** (see `AUDIT_LOG.md`, item 8); the values
above are from the corrected version.

## Not executed here

| script | why | what would be needed |
|---|---|---|
| `s01_build_ld.py` | needs the 1000G PLINK panel (~8.5 M variants) | the panel in `config.yaml` |
| `s02_extract_eqtl.py` | needs the bgzip/tabix GTEx files (4 GB per tissue) | the files listed under `eqtl` in `config.yaml` |
| `s03_dsi.py`, `s04_borrowing.py` | need the 10-trait MTAG output | `mtag_dir` |
| `s07`–`s09` | need GSCAN and the outcome GWAS | `gwas_dir` |
| `s11_mtag_validate.py` | needs a completed MTAG run | `mtag_dir`, `omega`/`sigma` |

These are the scripts that consume the raw inputs, which are public but large.
Everything they produce that the manuscript reports is shipped in
`source_data/`, and `verify_claims.py` re-derives the reported numbers from
those files, so the results can be checked without downloading them.
