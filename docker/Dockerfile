FROM quay.io/biocontainers/bioconductor-bambu:3.8.3--r44he5774e6_0
LABEL authors="preskaa"
# Install BiocManager from CRAN
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Optional: preinstall annotation packages used by bambu
# RUN Rscript -e "BiocManager::install('GenomicFeatures')"

# Optional: Set default command (if not using ENTRYPOINT for something else)
CMD ["R"]