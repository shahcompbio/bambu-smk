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
bambu_out <- snakemake@params[["outdir"]]
##prep annotations ...
fa.file <- snakemake@params[["ref_genome"]]
gtf.file <- snakemake@params[["ref_gtf"]]
yieldSize <- snakemake@params[["yieldsize"]]
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(gtf.file)
## load in all the bam files in current directory ...
## set the working directory
# outdir <- "/juno/work/shah/users/preskaa/APS012_Archive/"
# setwd(outdir)
#### check what RC outpaths exist and make a list
file_paths <- Sys.glob(file.path(bambu_out, "*", "rcOutput"))
print(file_paths)
###########
## make a vector of the read classes
library(BiocFileCache)
# Initialize an empty vector to contain the paths to the read classes ...
rds <- c()
# Loop through each file path
for (file_path in file_paths) {
  bfc <- BiocFileCache(file_path, ask = FALSE)
  info <- bfcinfo(bfc)
  ## grab last rpath in the cache; will figure out better system later
  # Append the rpath to the vector
  rds <- c(rds, info$rpath[length(info$rpath)])
}
## print read classes as well ...
print(rds)
####################
### run with recommended NDR
#############

## make output for multisample using recommended NDR
multi_out <- sprintf("%s/transcriptome_NDR_default/", bambu_out)
create_directory(multi_out)
ncores <- length(rds)
##trying to set baselineFDR vs NDR here to be more permissive ...
## baselineFDR defaults to 0.1 when NDR is not selected ...
## try with model training to begin with as Andre Sim recommended ...
se.NDR_default <- bambu(reads = rds, annotations = annotations, ncore = ncores,
            genome = fa.file, verbose=TRUE, yieldSize = yieldSize, lowMemory=FALSE,
            fusionMode=FALSE)

## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_default, path = multi_out)
#
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_default.RData", multi_out)
save(se.NDR_default, file = se_out)
### now let's test different NDR thresholds
## start by doing transcript discovery
newAnnotations <- bambu(reads = rds, annotations = annotations, ncore = ncores,
                        genome = fa.file, NDR = 1, verbose=TRUE,
                        yieldSize = yieldSize, lowMemory=FALSE, quant = FALSE)
##################
### do quantification on all new transcripts (NDR = 1) ....
###########
se.NDR_1 <- bambu(reads = rds, annotations = newAnnotations, ncore = ncores,
                  genome = fa.file, verbose=TRUE, yieldSize = yieldSize,
                  lowMemory=FALSE, NDR = 1, discovery = FALSE)
multi_out <- sprintf("%s/transcriptome_NDR_1/", bambu_out)
create_directory(multi_out)
## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_1, path = multi_out)
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_1.RData", multi_out)
save(se.NDR_1, file = se_out)
####################
## do quantification on NDR < 0.1 transcripts ...
##################
## make a filtered annotation
annotations.filtered <- newAnnotations[(!is.na(mcols(newAnnotations)$NDR) & mcols(newAnnotations)$NDR<0.1) | is.na(mcols(newAnnotations)$NDR)]
se.NDR_0.1 <- bambu(reads = rds, annotations = annotations.filtered, ncore = ncores,
                    genome = fa.file, verbose=TRUE, yieldSize = yieldSize,
                    lowMemory=FALSE, NDR = 1, discovery = FALSE)
multi_out <- sprintf("%s/transcriptome_NDR_0.1/", bambu_out)
create_directory(multi_out)
## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_0.1, path = multi_out)
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_0.1.RData", multi_out)
save(se.NDR_0.1, file = se_out)
