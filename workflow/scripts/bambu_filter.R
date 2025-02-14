#!/usr/bin/env Rscript

### define directory creation function
create_directory <- function(bambu_out) {
  if (!dir.exists(bambu_out)) {
    dir.create(bambu_out)
    print(paste("Directory created at", bambu_out))
  } else {
    print("Directory already exists")
  }
}

sepath <- snakemake@input[["sepath"]]
detectpath <- snakemake@output[["gtf"]]
n_samples <- snakemake@params[["n_samples"]]
##########
## prepare data
library(bambu)
##check version
# Assuming 'bambu' is already loaded
bambu_version <- packageVersion("bambu")
print(bambu_version)
## do work
load(sepath)
## rename the summarized experiment
se <- se.NDR_default
if (n_samples > 1) {
  fullLengthCounts <- assays(se)$fullLengthCounts
  rows_satisfying_conditions <- apply(fullLengthCounts, 1, function(row) any(row > 1))
  # Subset the SummarizedExperiment object based on the condition
  se.detected <- se[rows_satisfying_conditions, ]
} else {
  se.detected <- se[assays(se)$fullLengthCounts > 1,]
}
detectgtf <- rowRanges(se.detected)
writeToGTF(detectgtf, detectpath)