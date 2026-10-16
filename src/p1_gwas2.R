.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages({library(GenomicSEM); library(data.table)})
setwd("~/Cancer_GenomicSEM")
traits <- readLines("surviving_traits.txt")
loci <- fread("loci273.txt", header = FALSE)$V1
load("cov10_canonical.Rdata"); cov <- cov10
cat("traits:", paste(traits, collapse = ","), "\n")

cat("=== sumstats ===\n")
ss <- sumstats(files = paste0("mtag_in3/", traits, ".txt"),
               ref = "mini_ref_mtag.txt",
               trait.names = traits,
               se.logit = rep(TRUE, length(traits)),
               OLS = NULL,
               linprob = rep(FALSE, length(traits)),
               info.filter = 0.6, maf.filter = 0, keep.indel = FALSE,
               parallel = TRUE, cores = 10)
cat("sumstats rows:", nrow(ss), " cols:", ncol(ss), "\n")
cat("cols:", paste(colnames(ss), collapse = ","), "\n")
saveRDS(ss, "factor_sumstats.rds")

setDT(ss)
loci_present <- intersect(loci, ss$SNP)
cat("loci present in sumstats:", length(loci_present), "\n")
bg <- ss$SNP[seq(1, nrow(ss), by = 40)]
keep <- intersect(union(loci, bg), ss$SNP)
sub <- ss[SNP %in% keep]
cat("SNPs to model:", nrow(sub), "\n")

cat("=== commonfactorGWAS ===\n")
res <- commonfactorGWAS(covstruc = cov, SNPs = as.data.frame(sub),
                        estimation = "DWLS", parallel = TRUE, cores = 24,
                        smooth_check = FALSE, GC = "standard")
cat("cols:", paste(colnames(res), collapse = ","), "\n")
fwrite(res, "cf_gwas_loci.csv")
cat("WROTE cf_gwas_loci.csv rows", nrow(res), "\n")
cat("P1_GWAS_DONE\n")
