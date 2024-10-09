#!/usr/bin/env Rscript

setwd("./statGraph")

# Use system packages instead of RENV. Makes github workflow faster as it does not need to install R dependencies again.
if(nzchar(Sys.getenv("DISABLE_RENV"))) {
    renv::deactivate()
} else {
    renv::load()
}

devtools::test(stop_on_failure = TRUE)
devtools::build()
