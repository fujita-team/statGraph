# Use an appropriate base image
FROM ubuntu:22.04

# Set environment variables to non-interactive to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install dependencies
RUN apt-get update 
RUN apt install -y \
            r-base \
            libfreetype6-dev \
            libpng-dev \
            libtiff5-dev \
            libjpeg-dev \
            r-cran-tinytex \
            texlive \
            texinfo \
            texlive-fonts-extra \
            curl \
            libcurl4 \
            libharfbuzz-dev \
            libfribidi-dev \
            libcurl4-openssl-dev \
            libxml2-dev \
            libssl-dev \
            libfontconfig1-dev \
            libgit2-dev

# Update R packages
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('rARPACK')"
RUN Rscript -e "install.packages('foreach')"
RUN Rscript -e "install.packages('doParallel')"
RUN Rscript -e "install.packages('ks')"
RUN Rscript -e "install.packages('getopt')"
RUN Rscript -e "install.packages('renv')"
