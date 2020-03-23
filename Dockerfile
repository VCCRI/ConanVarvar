FROM rocker/shiny-verse:3.6.1

LABEL software = "ConanVarvar"
LABEL about.summary = "ConanVarvar: a versatile tool for the detection of large syndromic copy number variation from whole genome sequencing data"
LABEL container = "conanvarvar"

RUN apt-get update && apt-get install -my wget gnupg
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libbz2-dev \
    liblzma-dev \
    libicu-dev
RUN apt-get clean

RUN R -e "install.packages(c('BiocManager', 'optparse', 'foreach', 'doParallel', 'openxlsx', 'latex2exp', 'latticeExtra', 'shinythemes', 'shinyjs'))"
RUN R -e "BiocManager::install(c('Rsamtools', 'rCGH', 'CopywriteR', 'fastseg', 'exomeCopy'))"

COPY *.R /
COPY *.rds /
COPY *.md /

ARG GUI=false
ENV GUI=${GUI}

COPY conanvarvar.sh /
RUN chmod +x /conanvarvar.sh

EXPOSE 3838
ENTRYPOINT ./conanvarvar.sh $GUI $0 $@
