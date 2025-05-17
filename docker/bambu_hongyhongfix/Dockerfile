# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_21-R-4.5.0

# Install required Bioconductor package
RUN R -e 'BiocManager::install("BiocParallel")'
RUN R -e 'BiocManager::install("BiocFileCache")'
# Install remotes package to enable installing from GitHub
RUN R -e 'install.packages("remotes")'

# Install bambu from the specific GitHub branch
RUN R -e 'remotes::install_github("GoekeLab/bambu@HongYhong_fix")'
