# Use an appropriate base image
FROM ubuntu:22.04

# Set environment variables to non-interactive to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Copy the R script to the container
COPY install_packages.r /statgraph/install_packages.R

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
            libfontconfig1-dev

# Update R packages
RUN Rscript /statgraph/install_packages.R


