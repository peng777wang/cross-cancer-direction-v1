.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages({library(GenomicSEM); library(data.table)})
setwd("~/Cancer_GenomicSEM")
traits <- readLines("surviving_traits.txt")
loci <- fread("loci273.txt", header = FALSE)$V1
load("cov10_canonical.Rdata"); cov <- cov10

I <- cov$I
dimnames(I) <- list(traits, traits)
cat("== LDSC intercept matrix ==\n")
print(round(I, 3))
cat("min eigen(I):", min(eigen(I, symmetric = TRUE, only.values = TRUE)$values), "\n")

Z <- matrix(NA_real_, length(loci), length(traits))
for (i in seq_along(traits)) {
  d <- fread(sprintf("mtag_out/mtag_all_trait_%d.txt", i),
             select = c("SNP", "Z"))
  setkey(d, SNP); d <- d[!duplicated(SNP)]
  Z[, i] <- d[.(loci)]$Z
}
colnames(Z) <- traits

Iinv <- solve(I)
chi_gls <- apply(Z, 1, function(z) as.numeric(t(z) %*% Iinv %*% z))
p_gls <- pchisq(chi_gls, df = ncol(Z), lower.tail = FALSE)
chi_indep <- rowSums(Z^2)
p_indep <- pchisq(chi_indep, df = ncol(Z), lower.tail = FALSE)

# heterogeneity (Cochran Q on the z scale)
zbar <- rowMeans(Z)
Q <- rowSums((Z - zbar)^2)
Q_p <- pchisq(Q, df = ncol(Z) - 1, lower.tail = FALSE)

# discordant arm profile: traits whose sign opposes the majority sign
n_disc <- integer(length(loci))
med_disc <- numeric(length(loci))
med_conc <- numeric(length(loci))
for (i in seq_len(nrow(Z))) {
  z <- Z[i, ]
  maj <- if (sum(sign(z)) >= 0) 1 else -1
  disc <- sign(z) != maj
  n_disc[i] <- sum(disc)
  med_disc[i] <- if (any(disc)) median(abs(z[disc])) else NA_real_
  med_conc[i] <- if (any(!disc)) median(abs(z[!disc])) else NA_real_
}

dsi <- fread("dsi_loci.csv")
out <- data.table(SNP = loci, chi_gls = chi_gls, p_gls = p_gls,
                  chi_indep = chi_indep, p_indep = p_indep,
                  Q = Q, Q_p = Q_p, n_discordant = n_disc,
                  median_absz_discordant = med_disc,
                  median_absz_concordant = med_conc)
setkey(out, SNP); setkey(dsi, SNP)
m <- merge(dsi, out, by = "SNP")
m[, class := ifelse(DSI >= 0.999, "aligned",
                    ifelse(DSI <= -0.999, "reversed", "mixed"))]
fwrite(m, "p1_gls_loci.csv")

cat("\n== GLS omnibus by MTAG class ==\n")
md <- as.data.frame(m)
summ <- do.call(rbind, lapply(split(md, md$class), function(d) {
  data.frame(
    class = d$class[1],
    n = nrow(d),
    median_chi_gls = round(median(d$chi_gls), 1),
    frac_gls_p5e8 = round(mean(d$p_gls < 5e-8), 3),
    frac_gls_p05 = round(mean(d$p_gls < 0.05), 3),
    median_Q = round(median(d$Q), 1),
    median_n_discordant = as.numeric(median(d$n_discordant)),
    median_absz_disc = round(median(d$median_absz_discordant, na.rm = TRUE), 2),
    median_absz_conc = round(median(d$median_absz_concordant, na.rm = TRUE), 2),
    stringsAsFactors = FALSE)
}))
print(summ, row.names = FALSE)
fwrite(summ, "p1_gls_summary.csv")

cat("\n== reversed loci: discordant arm vs concordant arm ==\n")
rev <- m[class == "reversed"]
cat("n reversed:", nrow(rev), "\n")
cat("GLS p < 5e-8:", sum(rev$p_gls < 5e-8), "/", nrow(rev), "\n")
cat("discordant arm |z| < 1:", sum(rev$median_absz_discordant < 1, na.rm = TRUE),
    "/", nrow(rev), "\n")
cat("median |z| discordant:", round(median(rev$median_absz_discordant, na.rm = TRUE), 2),
    " concordant:", round(median(rev$median_absz_concordant, na.rm = TRUE), 2), "\n")
print(wilcox.test(rev$median_absz_discordant, rev$median_absz_concordant,
                  paired = TRUE))
cat("P1_GLS_DONE\n")
