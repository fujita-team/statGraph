#!/usr/bin/env Rscript
# This Rscript builds, test and constructs the documentation for statGraph
# It exists with 1 in care of error, and 0 in case of success
# It is invoked by GithubActions for automatic testing

if (!require("devtools")) install.packages("devtools")
if (!require("igraph")) install.packages("igraph")
if (!require("rARPACK")) install.packages("rARPACK")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")
if (!require("ks")) install.packages("ks")

library(devtools)

VALID_TESTS = list()
VALID_TESTS$mode <- "ALL"
PATH <- "./statGraph/"

# Run tests
devtools::test(pkg=PATH, stop_on_failure=TRUE)

# Generate Documentation
devtools::document(pkg=PATH)
devtools::build_manual(pkg=PATH)

# Build
devtools::build(pkg=PATH)
