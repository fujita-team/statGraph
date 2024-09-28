#!/usr/bin/Rscript
# This Rscript install all of the R dependencies for statGraph

if (!require("devtools")) install.packages("devtools")
if (!require("igraph")) install.packages("igraph")
if (!require("rARPACK")) install.packages("rARPACK")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")
if (!require("ks")) install.packages("ks")
if (!require("getopt")) install.packages("getopt")
