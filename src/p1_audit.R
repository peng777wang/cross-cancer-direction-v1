.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages({library(GenomicSEM); library(data.table)})
setwd("~/Cancer_GenomicSEM")
traits <- readLines("surviving_traits.txt")
loci <- fread("loci273.txt", header = FALSE)$V1
load("cov10_canonical.Rdata"); cov <- cov10

cat("=== (a) common factor model on the 10 traits (canonical cov) ===\n")
covs <- list(V = cov$V, S = cov$S, I = cov$I, N = cov$N, m = cov$m)
cf <- commonfactor(covstruc = covs, estimation = "DWLS")
print(cf$modelfit)
loadings <- cf$results[cf$results$op == "=~", ]
print(loadings[, c("lhs", "op", "rhs", "Unstandardized_Estimate",
                   "Standardized_Est", "p_value")])
fwrite(loadings, "p1_commonfactor_loadings.csv")

cat("\n=== (b) model-free multi-trait audit at the 273 loci ===\n")
z <- data.frame(SNP = loci, stringsAsFactors = FALSE)
n <- data.frame(SNP = loci, stringsAsFactors = FALSE)
for (i in seq_along(traits)) {
  d <- fread(sprintf("mtag_out/mtag_all_trait_%d.txt", i),
             select = c("SNP", "Z", "N"))
  setkey(d, SNP); d <- d[!duplicated(SNP)]
  r <- d[.(loci)]
  z[[traits[i]]] <- r$Z
  n[[traits[i]]] <- r$N
}
Z <- as.matrix(z[, traits])
N <- as.matrix(n[, traits])
w <- sqrt(N) / rowMeans(sqrt(N))

stouffer_equal <- rowSums(Z) / sqrt(ncol(Z))
stouffer_w <- rowSums(w * Z) / sqrt(rowSums(w^2))
zbar <- rowSums(w * Z) / rowSums(w)
Q <- rowSums(w * (Z - zbar)^2)
Q_p <- pchisq(Q, df = ncol(Z) - 1, lower.tail = FALSE)
p_two <- 2 * pnorm(-abs(Z))
minP <- apply(p_two, 1, min)
n_sig <- rowSums(p_two < 5e-8)
n_nom <- rowSums(p_two < 0.05)

audit <- data.table(SNP = loci, stouffer_equal = stouffer_equal,
                    stouffer_w = stouffer_w, Q = Q, Q_p = Q_p,
                    minP = minP, n_sig = n_sig, n_nom = n_nom)
dsi <- fread("dsi_loci.csv")
setkey(audit, SNP); setkey(dsi, SNP)
m <- merge(dsi, audit, by = "SNP", all.x = TRUE)
m[, class := ifelse(DSI >= 0.999, "aligned",
                    ifelse(DSI <= -0.999, "reversed", "mixed"))]
fwrite(m, "p1_locus_audit.csv")
cat("classes:\n"); print(m[, .N, by = class])

cat("\n=== (c) by MTAG direction class ===\n")
summ <- m[, .(n = .N,
              median_abs_stouffer = round(median(abs(stouffer_w)), 2),
              frac_stouffer_p5e8 = round(mean(2 * pnorm(-abs(stouffer_w)) < 5e-8), 3),
              median_Q = round(median(Q), 1),
              frac_Q_p05 = round(mean(Q_p < 0.05), 3),
              median_minP = signif(median(minP), 3),
              median_n_sig = median(n_sig)), by = class]
print(summ)
fwrite(summ, "p1_class_summary.csv")
cat("P1_AUDIT_DONE\n")
