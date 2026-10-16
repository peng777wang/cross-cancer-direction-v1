suppressPackageStartupMessages(library(data.table))
setwd("~/Cancer_GenomicSEM")
res <- fread("pleiotropy_resource.csv")
cat("rows", nrow(res), "loci", uniqueN(res$snp), "traits",
    uniqueN(res$trait), "\n")
cat("classes:\n"); print(res[, .N, by = class_mtag])
rev <- res[grepl("reversed", class_mtag, fixed = TRUE)]
cat("\n== reversed class: MTAG-significant rows by trait ==\n")
print(rev[mtag_significant == TRUE,
          .(n = .N,
            median_borrowed = round(median(borrowed_fraction), 4),
            median_abs_zraw = round(median(abs(z_raw)), 3),
            own_gws = sum(own_gwas_significant)),
          by = trait][order(-n)])
cat("\n== reversed loci: k distribution ==\n")
print(rev[, .(k = sum(mtag_significant)), by = snp][, .N, by = k])
cat("\n== aligned loci: k distribution ==\n")
al <- res[class_mtag == "aligned"]
print(al[, .(k = sum(mtag_significant)), by = snp][, .N, by = k])
cat("\n== PC rows inside the reversed class (any significance flag) ==\n")
print(rev[trait == "PC", .(n = .N, n_sig = sum(mtag_significant),
                           median_borrowed = round(median(borrowed_fraction), 4))])
cat("\n== reversed class: own-GWAS significant cancers per locus ==\n")
print(rev[, .(own = sum(own_gwas_significant)), by = snp][, .N, by = own])
cat("\nP1_REV_CHECK_DONE\n")
