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
if (!require("getopt")) install.packages("getopt")


library(devtools)
library(getopt)

# Parse command line arguments

get_list <- function(string){
    output <- strsplit(string, ',')[[1]]
    return(output)
}
opt_all <- c("all", "a", 0, "character", "Perform all tests")
opt_include <- c("include", "i", 1, "character", "Comma separated list of tests to include")
opt_exclude <- c("exclude", "e", 1, "character", "Comma separated list of tests to exclude")
opt_help <- c("help", "h", 0, "character", "Prints the help menu.")

spec <- rbind(opt_all, opt_include, opt_exclude, opt_help)
opt <- getopt(spec)

if(!is.null(opt$help)){
    
    cat("Automatically performs tests for the statGraph package.\n\n")
    cat("Usage: ./test.r [ARGUMENTS]\n\n")
    cat("The following is the valid set of possible arguments\n")
    cat("-a\t\t\t\t-\tRuns all tests\n")
    cat("-i TEST1,TEST2,TEST2\t\t-\tA comma separated list of test names including the names of the tests that should be included\n")
    cat("-e TEST1,TEST2,TEST2\t\t-\tA comma separated list of test names including the names of the tests that should be excluded\n")
    cat("-h\t\t\t\t-\tPrints this menu\n")
    quit("yes")
}

VALID_TESTS <- list()

if(is.null(opt$all) && is.null(opt$include) && is.null(opt$exclude)){
    VALID_TESTS$mode <- "ALL"
} else {
    if(!is.null(opt$all)){
        VALID_TESTS$mode <- "ALL"
    } 
    else if(!is.null(opt$include)) {
        VALID_TESTS$mode <- "WHITELIST"
        VALID_TESTS$list <- get_list(opt$include)
    }
    else if(!is.null(opt$exclude)) {
        VALID_TESTS$mode <- "BLACKLIST"
        VALID_TESTS$list <- get_list(opt$exclude)
    }
}

# Run Tests

PATH <- "./statGraph/"

# Run tests
devtools::test(pkg=PATH, stop_on_failure=FALSE)

