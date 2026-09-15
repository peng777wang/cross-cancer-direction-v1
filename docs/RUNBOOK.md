# Runbook: the exact command for each step

Commands are written against the paths declared in `config.yaml`. Nothing here
touches a server; the scripts are self-contained and read their inputs from the
paths you pass.

Steps marked **[executed]** were run and their output recorded in
`EXECUTION_LOG.md`. Steps marked **[not executed here]** consume the large public
inputs and were not re-run in this environment; their outputs are shipped in
`source_data/` so that the reported numbers remain checkable.

## Step 0 — what you can run immediately

```bash
python src/verify_claims.py                 # [executed] 19/19 PASS
python tests/test_reproducibility.py        # [executed] 22 tests OK
python src/s05_coloc_core.py                # [executed] four-scenario self-test
python src/s12_mtag_simulation.py --reps 20 --snps 200000 --out outputs
```

## Step 1 — LD matrix for a region

```bash
python src/s01_build_ld.py \
  --panel /path/EUR_local/EUR --chrom 15 --start 78400000 --end 79300000 \
  --gwas-dir /path/clean --traits LUAD LUSC SCLC --out outputs
```

## Step 2 — eQTL z-scores for the region, harmonised to the reference alleles

```bash
python src/s02_extract_eqtl.py \
  --eqtl /path/Lung.tsv.gz --chrom 15 --start 78600000 --end 79100000 \
  --gene ENSG00000136381 IREB2 --meta outputs/ld_meta.tsv --out outputs
```

## Step 3 — direction architecture

```bash
python src/s03_dsi.py --mtag-dir /path/mtag_out --out outputs          # [not executed here]
python src/s04_borrowing.py --mtag-dir /path/mtag_out --out outputs    # [not executed here]
```

These two are the retrospective half of the paper: `s03` reproduces the aligned
and reversed classes (and both of their null comparisons), `s04` shows how much
of each reversed-arm estimate is borrowed rather than observed.

## Step 4 — colocalisation

```bash
Rscript src/s06_coloc_susie.R \
  outputs/ld.rds outputs/ld_meta.tsv outputs/gwas_z.tsv \
  outputs/eqtl_z.tsv outputs/n.txt outputs/coloc_out.csv
```

**[executed]** on the cluster: IREB2 x signal B gave PP.H4 = 0.715 / 0.962 /
0.996 at p12 = 1e-6 / 1e-5 / 1e-4, and IREB2 x signal A gave PP.H3 >= 0.82.

Two things must stay true for this script to be valid:

* the LD matrix's row order must equal the metadata table's SNP order; the script
  asserts this (`stopifnot(all(rownames(R) == meta$snp))`);
* the LD matrix must be indexed with positions in the **full** metadata table —
  subsetting the table first and then indexing the matrix silently mismatches.

## Step 5 — instruments and MR

```bash
python src/s07_instruments.py --gscan /path/GSCAN.txt.gz --panel /path/EUR \
  --chrom 15 11 2 1 --pos 78600319 113040282 145385522 43572014

python src/s08_cis_mr.py --beta-eqtl -0.108214 --se-eqtl 0.0203518 \
  --outcome LUAD -0.225095 0.0194132

python src/s09_mvmr.py      # exposes mvmr(), wald_ratio(), marginal_f(), cochran_q()
```

## Step 6 — MTAG port and its validation

```bash
python src/s11_mtag_validate.py \
  --omega /path/mtag_all_omega_hat.txt --sigma /path/mtag_all_sigma_hat.txt \
  --mtag-dir /path/mtag_out --traits LUAD LUSC SCLC EC EAC THCA RCC PC GC PRAD \
  --snps 20000 --out outputs
```

Expected: Pearson r = 1.00000000 against the official MTAG output and a maximum
absolute difference of order 1e-15.
