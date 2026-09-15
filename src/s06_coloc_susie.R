# susie_rss + coloc.susie at a single region, over several p12 priors.
#
# Usage: Rscript s06_coloc_susie.R <ld_rds> <meta_tsv> <gwas_tsv> <eqtl_tsv> <n_file> <out_csv>
#   ld_rds    : correlation matrix (RDS), rownames/colnames = rsIDs
#   meta_tsv  : snp, a1, a2, bp, f_ref_A1          (rows define the SNP universe)
#   gwas_tsv  : snp, z_gwas                        (harmonised to meta a1)
#   eqtl_tsv  : snp, <gene1>, <gene2>, ...         (harmonised to meta a1)
#   n_file    : two lines: n_gwas, n_eqtl
#
# IMPORTANT: the LD matrix must be indexed with positions in the FULL metadata
# table. Subsetting the table first and then indexing the matrix is a silent
# mismatch (see docs/AUDIT_LOG.md item 7).

suppressMessages({library(susieR); library(coloc)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("six arguments required")
ld_rds <- args[1]; meta_f <- args[2]; gwas_f <- args[3]
eqtl_f <- args[4]; n_f <- args[5]; out_f <- args[6]

WIN <- c(78800000, 79000000)          # hg19 window; edit for other loci
P12 <- c(1e-6, 1e-5, 1e-4)

R    <- readRDS(ld_rds)
meta <- read.table(meta_f, header = TRUE, stringsAsFactors = FALSE)
gw   <- read.table(gwas_f, header = TRUE, stringsAsFactors = FALSE)
eq   <- read.table(eqtl_f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
nn   <- scan(n_f, quiet = TRUE)

stopifnot(all(rownames(R) == meta$snp))      # the matrix must follow the table
snps <- meta$snp[meta$bp >= WIN[1] & meta$bp <= WIN[2]]

res <- list()
for (gene in setdiff(colnames(eq), "snp")) {
  d <- data.frame(snp = snps)
  d$z_gwas <- gw$z_gwas[match(snps, gw$snp)]
  d$z_eqtl <- eq[[gene]][match(snps, eq$snp)]
  d <- d[!is.na(d$z_gwas) & !is.na(d$z_eqtl), ]
  idx <- match(d$snp, meta$snp)               # indices into the FULL matrix
  Rs <- R[idx, idx]
  Rs <- (Rs + t(Rs)) / 2
  ev <- min(eigen(Rs, symmetric = TRUE, only.values = TRUE)$values)
  if (ev < 1e-4) Rs <- Rs + diag(1e-4 - ev, nrow(Rs))

  s1 <- susie_rss(z = d$z_eqtl, R = Rs, n = nn[2], L = 5, max_iter = 1000)
  s2 <- susie_rss(z = d$z_gwas, R = Rs, n = nn[1], L = 5, max_iter = 1000)
  if (!isTRUE(s1$converged) || !isTRUE(s2$converged)) {
    message("  ", gene, ": susie did not converge; skipped")
    next
  }
  for (p12 in P12) {
    cc <- try(coloc.susie(s1, s2, p12 = p12), silent = TRUE)
    # coloc.susie returns list(nsnps = ...) -- with NO $summary element -- when
    # neither trait contributes a credible set, or when no credible-set pair
    # passes the p12 threshold.  Guard on the element, not on the object.
    if (inherits(cc, "try-error") || is.null(cc) || is.null(cc$summary)) next
    sm <- as.data.frame(cc$summary, stringsAsFactors = FALSE)
    if (nrow(sm) == 0) next
    sm$gene <- gene; sm$p12 <- p12; sm$n_snp <- nrow(d)
    res[[paste(gene, p12)]] <- sm
  }
}
if (length(res)) {
  out <- do.call(rbind, res)
  write.csv(out, out_f, row.names = FALSE)
  keep <- intersect(c("gene", "p12", "n_snp", "idx1", "idx2", "hit1", "hit2",
                      "PP.H3.abf", "PP.H4.abf"), colnames(out))
  print(out[, keep, drop = FALSE])
} else {
  message("no colocalised pair in any gene / prior")
}
