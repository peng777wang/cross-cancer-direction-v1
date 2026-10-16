.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("~/Cancer_GenomicSEM")
traits <- readLines("surviving_traits.txt")
cat("traits:", length(traits), "\n")
z <- NULL
for (i in seq_along(traits)) {
  d <- fread(sprintf("mtag_out/mtag_all_trait_%d.txt", i),
             select = c("SNP", "Z"))
  setkey(d, SNP)
  d <- d[!duplicated(SNP)]
  setnames(d, "Z", traits[i])
  z <- if (is.null(z)) d else merge(z, d, by = "SNP")
  cat("  ", i, traits[i], nrow(z), "\n")
}
m <- as.matrix(z[, ..traits])
cat("SNPs:", nrow(m), "\n")
R_all <- cor(m, use = "pairwise.complete.obs")
keep <- rowSums(abs(m) > 5.451) == 0
cat("SNPs with no trait above 5.451:", sum(keep), "\n")
R_null <- cor(m[keep, ], use = "pairwise.complete.obs")
write.csv(round(R_all, 4), "p1_trait_cor_all.csv")
write.csv(round(R_null, 4), "p1_trait_cor_null.csv")
cat("\n== correlation among the ten MTAG z-scores, all SNPs ==\n")
print(round(R_all, 3))
cat("\n== off-diagonal summary, all SNPs: min", round(min(R_all[upper.tri(R_all)]), 3),
    " median", round(median(R_all[upper.tri(R_all)]), 3),
    " max", round(max(R_all[upper.tri(R_all)]), 3), "\n")
cat("== off-diagonal summary, null SNPs: min", round(min(R_null[upper.tri(R_null)]), 3),
    " median", round(median(R_null[upper.tri(R_null)]), 3),
    " max", round(max(R_null[upper.tri(R_null)]), 3), "\n")
cat("P1_TRAIT_COR_DONE\n")
