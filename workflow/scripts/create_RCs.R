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

##########
## prepare data
library(bambu)
##check version
# Assuming 'bambu' is already loaded
bambu_version <- packageVersion("bambu")
print(bambu_version)

## load in commandline args
sample_name <- snakemake@params[["sample_name"]]
sample.bam <- snakemake@input[["bam"]]
bambu_out <- snakemake@params[["outdir"]]
yieldSize <- snakemake@params[["yieldsize"]]
##prep annotations ...
fa.file <- snakemake@params[["ref_genome"]]
gtf.file <- snakemake@params[["ref_gtf"]]
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(gtf.file)
## load in all the bam files in current directory ...
## set the working directory

## make bambu outdir
create_directory(bambu_out)
##create output directory for rcFiles
rcOut <- sprintf("%s/%s/rcOutput/", bambu_out, sample_name)
create_directory(rcOut)
## delete any screwed up previous RCfiles

# List all files in the directory with the ".rds" extension
rds_files <- list.files(rcOut, pattern = "\\.rds$", full.names = TRUE)

# Delete the files
if (length(rds_files) > 0) {
  unlink(rds_files)
  cat("Files with .rds extension deleted successfully.\n")
} else {
  cat("No files with .rds extension found in the directory.\n")
}

## create RC files
se.discoveryOnly <- bambu(reads = sample.bam, annotations = annotations, ncore = 1, genome = fa.file,
            yieldSize = yieldSize, lowMemory=TRUE, fusionMode=FALSE,
            verbose=TRUE, discovery=TRUE, quant=FALSE, rcOutDir = rcOut)
## make an extended gtf for the sample, we will do quantification in following rule; this is mostly to give snakemake
## a file to target
outgtf <- sprintf("%s/%s/output.gtf", bambu_out, sample_name)
writeToGTF(se.discoveryOnly, outgtf)
