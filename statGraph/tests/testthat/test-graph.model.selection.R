# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

# RunTest("graph.model.selection", {
#             G <- igraph::sample_gnp(n=30, p=0.5)
#             #' # Erdos-Renyi graph
#             model1 <- function(n, p){
#                 return(igraph::sample_gnp(n, p))
#             }
#             # Watts-Strogatz small-world graph
#             model2 <- function(n, pr, K=8){
#                 return(igraph::sample_smallworld(1, n, K, pr))
#             }
#             parameters <- list(seq(0.1, 0.9, 0.1), seq(0.1, 0.9, 0.1))
#             result2 <- graph.model.selection(G, list(model1, model2), parameters)
#             result2
#
#
#             test_that("graph.model.selection for undirected graphs (ternary search)", {
#                           expect_equal(result2$model, 1)
# })
# })
#
# RunTest("graph.model.selection", {
#             #' # Erdos-Renyi graph
#             model1 <- function(n, p){
#                 U <- as.matrix( igraph::as_adj( igraph::sample_gnp(n, p)))
#                 L <- as.matrix( igraph::as_adj( igraph::sample_gnp(n, p)))
#                 U[upper.tri(U)] <- 0
#                 L[lower.tri(L)] <- 0
#
#                 return(U + L)
#             }
#             # Watts-Strogatz small-world graph
#             model2 <- function(n, p, K=8){
#                 U <- as.matrix( igraph::as_adj( igraph::sample_smallworld(1, n, K, p)))
#                 L <- as.matrix( igraph::as_adj( igraph::sample_smallworld(1, n, K, p)))
#                 U[upper.tri(U)] <- 0
#                 L[lower.tri(L)] <- 0
#
#                 return(U + L)
#             }
#             G <- model1(30, 0.3)
#
#             parameters <- list(seq(0.1, 0.9, 0.1), seq(0.1, 0.9, 0.1))
#             result2 <- graph.model.selection(G, list(model1, model2), parameters, directed=TRUE)
#             result2
#
#
#             test_that("graph.model.selection for directed graphs (grid search)", {
#                           expect_equal(result2$model, 1)
# })
# })
