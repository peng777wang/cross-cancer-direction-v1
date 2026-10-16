.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("~/Cancer_GenomicSEM")
res <- fread("pleiotropy_resource.csv")
cat("rows:", nrow(res), " loci:", uniqueN(res$snp), "\n")
cat("class_mtag:\n"); print(table(res$class_mtag))
sig <- res[mtag_significant == TRUE]
sig[, arm := ifelse(sign(z_mtag) != sign(z_raw), "discordant", "concordant")]
cat("\nMTAG-significant rows:", nrow(sig), "\n")
print(table(sig$class_mtag, sig$arm))
cat("\n== among MTAG-significant rows, by arm ==\n")
print(sig[, .(n = .N,
              median_abs_zraw = round(median(abs(z_raw)), 3),
              median_abs_zmtag = round(median(abs(z_mtag)), 3),
              median_borrowed = round(median(borrowed_fraction), 4),
              frac_own_gws = round(mean(own_gwas_significant), 4)),
          by = .(class_mtag, arm)][order(class_mtag, arm)])
cat("\n== paired comparison at the reversed loci ==\n")
rev <- sig[grepl("reversed", class_mtag, fixed = TRUE)]
locs <- split(rev, rev$snp)
pair <- rbindlist(lapply(locs, function(d) {
  dis <- d[arm == "discordant"]
  con <- d[arm == "concordant"]
  if (nrow(dis) == 0 || nrow(con) == 0) return(NULL)
  data.table(snp = d$snp[1], n_disc = nrow(dis), n_conc = nrow(con),
             max_abs_zraw_disc = max(abs(dis$z_raw)),
             max_abs_zraw_conc = max(abs(con$z_raw)),
             max_borrowed_disc = max(dis$borrowed_fraction),
             any_own_gws_conc = any(con$own_gwas_significant),
             any_own_gws_disc = any(dis$own_gwas_significant))
}))
cat("reversed loci with both arms:", nrow(pair), "\n")
if (nrow(pair) > 0) {
  cat("discordant arm own GWS:", sum(pair$any_own_gws_disc), "of", nrow(pair), "\n")
  cat("concordant arm own GWS:", sum(pair$any_own_gws_conc), "of", nrow(pair), "\n")
  cat("median max|z_raw| discordant:", round(median(pair$max_abs_zraw_disc), 3),
      " concordant:", round(median(pair$max_abs_zraw_conc), 3), "\n")
  print(wilcox.test(pair$max_abs_zraw_disc, pair$max_abs_zraw_conc, paired = TRUE))
  cat("median max borrowed (discordant arm):",
      round(median(pair$max_borrowed_disc), 4), "\n")
  fwrite(pair, "p1_reversed_arm_pairs.csv")
}
cat("\n== control: aligned loci ==\n")
aligned_rows <- sig[class_mtag == "aligned"]
cat("aligned significant rows:", nrow(aligned_rows), " discordant:",
    sum(aligned_rows$arm == "discordant"), "\n")
cat("aligned own GWS rate:", round(mean(aligned_rows$own_gwas_significant), 3), "\n")
cat("P1_ARMS_DONE\n")
