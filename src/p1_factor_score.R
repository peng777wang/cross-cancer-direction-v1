.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages({library(GenomicSEM); library(data.table)})
setwd("~/Cancer_GenomicSEM")
traits <- readLines("surviving_traits.txt")
loci <- fread("loci273.txt", header = FALSE)$V1
load("cov10_canonical.Rdata"); cov <- cov10

S <- cov$S
dimnames(S) <- list(traits, traits)
eg <- eigen(S)
w1 <- eg$vectors[, 1]
if (sum(w1) < 0) w1 <- -w1
cat("leading eigenvalue share:", round(eg$values[1] / sum(eg$values), 3), "\n")
cat("factor weights:\n")
print(round(setNames(w1, traits), 3))

Z <- matrix(NA_real_, length(loci), length(traits))
N <- matrix(NA_real_, length(loci), length(traits))
for (i in seq_along(traits)) {
  d <- fread(sprintf("mtag_out/mtag_all_trait_%d.txt", i),
             select = c("SNP", "Z", "N"))
  setkey(d, SNP); d <- d[!duplicated(SNP)]
  r <- d[.(loci)]
  Z[, i] <- r$Z; N[, i] <- r$N
}
Nfac <- N / rowMeans(N, na.rm = TRUE)
denom <- sqrt(rowSums((w1^2) * Nfac))
score <- as.vector(Z %*% w1) / denom

dsi <- fread("dsi_loci.csv")
out <- data.table(SNP = loci, factor_score_z = score,
                  factor_score_p = 2 * pnorm(-abs(score)))
setkey(out, SNP); setkey(dsi, SNP)
out <- merge(out, dsi, by = "SNP")
out[, class := ifelse(DSI >= 0.999, "aligned",
                      ifelse(DSI <= -0.999, "reversed", "mixed"))]
fwrite(out, "p1_factor_score.csv")
cat("\nby class:\n")
print(out[, .(n = .N,
              median_score = round(median(abs(factor_score_z)), 2),
              frac_p5e8 = round(mean(factor_score_p < 5e-8), 3),
              frac_p01 = round(mean(factor_score_p < 0.05), 3)), by = class])
cat("P1_FACTOR_DONE\n")
