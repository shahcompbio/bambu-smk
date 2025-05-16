# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_21-R-4.5.0

# Install required Bioconductor package
RUN R -e 'BiocManager::install("BiocParallel")'
RUN R -e 'BiocManager::install("BiocFileCache")'
RUN R -e 'BiocManager::install("bambu")'