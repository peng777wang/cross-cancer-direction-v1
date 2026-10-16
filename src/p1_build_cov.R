.libPaths(c(path.expand("~/Rlib_gsem"), .libPaths()))
suppressPackageStartupMessages(library(GenomicSEM))
setwd("~/Cancer_GenomicSEM")
surv <- readLines("surviving_traits.txt")
pop.prev <- c(LUAD = 0.040, LUSC = 0.020, SCLC = 0.005, EC = 0.030,
              EAC = 0.005, THCA = 0.012, RCC = 0.018, PC = 0.016,
              GC = 0.006, PRAD = 0.127)
cat("traits:", paste(surv, collapse = ","), "\n")
cat("population prevalences:", paste(pop.prev[surv], collapse = ","), "\n")
cov10 <- ldsc(traits = paste0(surv, ".sumstats.gz"),
              sample.prev = rep(0.5, length(surv)),
              population.prev = pop.prev[surv],
              ld = "ldref/eur_w_ld_chr/",
              wld = "ldref/eur_w_ld_chr/",
              trait.names = surv)
save(cov10, file = "cov10_canonical.Rdata")
cat("S:", paste(dim(cov10$S), collapse = "x"),
    " V:", paste(dim(cov10$V), collapse = "x"),
    " I:", paste(dim(cov10$I), collapse = "x"),
    " N:", paste(dim(cov10$N), collapse = "x"), "\n")
cat("h2:\n"); print(cov10$S)
cat("COV_BUILD_DONE\n")
