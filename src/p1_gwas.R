.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})
cat("GenomicSEM from:", find.package("GenomicSEM"), "\n")
setwd("~/Cancer_GenomicSEM")
args <- commandArgs(trailingOnly = TRUE)
TEST <- length(args) > 0 && args[1] == "test"
NCORE <- 24

traits <- readLines("surviving_traits.txt")
loci <- fread("loci273.txt", header = FALSE)$V1
cat("traits:", paste(traits, collapse = ","), "\n")

COVFILE <- if (file.exists("cov10_canonical.Rdata")) "cov10_canonical.Rdata" else "ldsc.cov.10.resourced.Rdata"
cat("using covariance file:", COVFILE, "\n")
load(COVFILE)
if (exists("cov10")) cov <- cov10
cat("cov names:", paste(names(cov), collapse = ","), "\n")
cat("S:", paste(dim(cov$S), collapse = "x"), " V:", paste(dim(cov$V), collapse = "x"),
    " I:", paste(dim(cov$I), collapse = "x"), " N:", paste(dim(cov$N), collapse = "x"), "\n")

d0 <- fread("mtag_out/mtag_all_trait_1.txt",
            select = c("SNP", "A1", "A2", "Z", "N", "FRQ"))
cat("MTAG output SNPs:", nrow(d0), "\n")
cat("loci present in MTAG output:", sum(loci %in% d0$SNP), "\n")
setkey(d0, SNP)
d0 <- d0[!duplicated(SNP)]
bg <- d0$SNP[seq(1, nrow(d0), by = 40)]
keep <- intersect(union(loci, bg), d0$SNP)
if (TEST) keep <- head(keep, 60)
cat("SNPs to model:", length(keep), "\n")

base <- d0[.(keep)]
snps <- data.frame(SNP = base$SNP, A1 = base$A1, A2 = base$A2,
                   MAF = base$FRQ, stringsAsFactors = FALSE)
flips <- integer(length(traits))
for (i in seq_along(traits)) {
  tr <- traits[i]
  d <- fread(sprintf("mtag_out/mtag_all_trait_%d.txt", i),
             select = c("SNP", "A1", "A2", "Z", "N", "FRQ"))
  setkey(d, SNP)
  d <- d[!duplicated(SNP)]
  d <- d[.(keep)]
  flipped <- (d$A1 == base$A2) & (d$A2 == base$A1)
  flips[i] <- sum(flipped)
  zz <- d$Z
  zz[flipped] <- -zz[flipped]
  varSNP <- 2 * d$FRQ * (1 - d$FRQ)
  varSNP[varSNP <= 0] <- NA
  denom <- sqrt((d$N / 4) * varSNP)
  snps[[paste0("beta.", tr)]] <- zz / denom
  snps[[paste0("se.", tr)]] <- 1 / denom
  snps[[paste0("N.", tr)]] <- d$N
}
cat("allele flips per trait:", paste(flips, collapse = ","), "\n")
chk <- abs(snps[["beta.LUAD"]] / snps[["se.LUAD"]] - base$Z)
cat("max |beta/se - Z| (LUAD):", max(chk, na.rm = TRUE), "\n")
cat("NA betas:", sum(is.na(snps[["beta.LUAD"]])), "\n")

cat("running commonfactorGWAS on", nrow(snps), "SNPs\n")
res <- commonfactorGWAS(covstruc = cov, SNPs = snps, estimation = "DWLS",
                        parallel = TRUE, cores = NCORE, smooth_check = FALSE,
                        GC = "standard")
cat("columns:", paste(colnames(res), collapse = ","), "\n")
out <- if (TEST) "cf_gwas_test.csv" else "cf_gwas_loci.csv"
fwrite(res, out)
cat("WROTE", out, "rows", nrow(res), "\n")
