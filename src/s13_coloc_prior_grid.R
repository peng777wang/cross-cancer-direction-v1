# LIMITATION 1: "PP.H4 depends on the p12 prior."
#
# The H3/H4 decision depends on the priors only through the single ratio
#     r = p12 / (p1 * p2)
# because PP.H4 ~ p12 * sum(e1*e2) and PP.H3 ~ p1*p2 * (sum e1 * sum e2 - sum e1*e2).
# So instead of quoting three numbers, this script maps PP.H4 over the whole
# plausible prior grid and reports the exact tipping point of r.
#
# Reuses the susie fits once per gene and re-evaluates coloc.susie for every
# (p1=p2, p12) combination, so the credible sets are held fixed and only the
# priors vary.

# If coloc/susieR are installed in a personal library, set R_LIBS_EXTRA.
if (nzchar(Sys.getenv("R_LIBS_EXTRA"))) .libPaths(c(Sys.getenv("R_LIBS_EXTRA"), .libPaths()))
suppressMessages({library(susieR); library(coloc)})

ROOT <- Sys.getenv("CGS_ROOT", unset = "data")   # path to the region inputs
OUT <- file.path(ROOT, "lim1_out")
dir.create(OUT, showWarnings = FALSE)

WIN <- c(78800000, 79000000)
P1P2 <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
P12 <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4)

R <- readRDS(file.path(ROOT, "release_ld.rds"))
meta <- read.table(file.path(ROOT, "release_meta.tsv"), header = TRUE,
                   stringsAsFactors = FALSE)
gw <- read.table(file.path(ROOT, "release_gwas.tsv"), header = TRUE,
                 stringsAsFactors = FALSE)
eq <- read.table(file.path(ROOT, "release_eqtl.tsv"), header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
nn <- scan(file.path(ROOT, "release_n.txt"), quiet = TRUE)
stopifnot(all(rownames(R) == meta$snp))

snps <- meta$snp[meta$bp >= WIN[1] & meta$bp <= WIN[2]]
rows <- list()

for (gene in setdiff(colnames(eq), "snp")) {
  d <- data.frame(snp = snps)
  d$z_gwas <- gw$z_gwas[match(snps, gw$snp)]
  d$z_eqtl <- eq[[gene]][match(snps, eq$snp)]
  d <- d[!is.na(d$z_gwas) & !is.na(d$z_eqtl), ]
  idx <- match(d$snp, meta$snp)
  Rs <- R[idx, idx]
  Rs <- (Rs + t(Rs)) / 2
  ev <- min(eigen(Rs, symmetric = TRUE, only.values = TRUE)$values)
  if (ev < 1e-4) Rs <- Rs + diag(1e-4 - ev, nrow(Rs))

  s1 <- susie_rss(z = d$z_eqtl, R = Rs, n = nn[2], L = 5, max_iter = 1000)
  s2 <- susie_rss(z = d$z_gwas, R = Rs, n = nn[1], L = 5, max_iter = 1000)
  if (!isTRUE(s1$converged) || !isTRUE(s2$converged)) next

  for (p in P1P2) {
    for (q in P12) {
      cc <- try(coloc.susie(s1, s2, p1 = p, p2 = p, p12 = q), silent = TRUE)
      if (inherits(cc, "try-error") || is.null(cc$summary)) next
      sm <- as.data.frame(cc$summary)
      if (nrow(sm) == 0) next
      sm$gene <- gene
      sm$p1 <- p
      sm$p12 <- q
      sm$prior_ratio <- q / (p * p)
      sm$n_snp <- nrow(d)
      rows[[length(rows) + 1]] <- sm
    }
  }
  cat("done", gene, "\n")
}

if (length(rows)) {
  out <- do.call(rbind, rows)
  keep <- intersect(c("gene", "p1", "p12", "prior_ratio", "idx1", "idx2",
                      "hit1", "hit2", "PP.H3.abf", "PP.H4.abf", "n_snp"),
                    colnames(out))
  out <- out[, keep]
  write.csv(out, file.path(OUT, "lim1_prior_grid.csv"), row.names = FALSE)
  cat("rows:", nrow(out), "\n")
} else {
  cat("no credible-set pairs\n")
}
