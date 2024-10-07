# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

RunTest("graph.model.selection", {
            #' # Erdos-Renyi graph
            model1 <- function(n, p){
                return(igraph::sample_gnp(n, p))
            }
            # Watts-Strogatz small-world graph
            model2 <- function(n, pr, K=8){
                return(igraph::sample_smallworld(1, n, K, pr))
            }

            G <- model1(n=30, p=0.5)
            parameters <- list(seq(0.1, 0.9, 0.1), seq(0.1, 0.9, 0.1))
            result2 <- graph.model.selection(G, list(model1, model2), parameters)

            choice <- as.numeric(which.min(result2$estimate[, 2]))
            test_that("graph.model.selection is working for undirected graphs", {
                          expect_equal(choice, 1)
            })
})

RunTest("graph.model.selection", {
            #' # Erdos-Renyi graph
            model1 <- function(n, p){
                return(igraph::sample_gnp(n, p))
            }
            # Watts-Strogatz small-world graph
            model2 <- function(n, pr, K=8){
                return(igraph::sample_smallworld(1, n, K, pr))
            }

            G <- model2(n=30, p=0.5)
            parameters <- list(seq(0.1, 0.9, 0.1), seq(0.1, 0.9, 0.1))
            result2 <- graph.model.selection(G, list(model1, model2), parameters)

            choice <- as.numeric(which.min(result2$estimate[, 2]))
            test_that("graph.model.selection is working for undirected graphs", {
                          expect_equal(choice, 2)
            })
})


RunTest("graph.model.selection", {
            #' # Erdos-Renyi graph
            model1 <- function(n, p){
                return(igraph::sample_gnp(n, p))
            }
            # Watts-Strogatz small-world graph
            model2 <- function(n, pr, K=8){
                return(igraph::sample_smallworld(1, n, K, pr))
            }

            G <- model1(n=30, p=0.5)
            parameters <- list(seq(0.1, 0.9, 0.01), seq(0.1, 0.9, 0.01))
            result2 <- graph.model.selection(G, list(model1, model2), parameters)

            choice <- as.numeric(which.min(result2$estimate[, 2]))
            test_that("graph.model.selection is working for undirected graphs", {
                          expect_equal(choice, 1)
            })
})

RunTest("graph.model.selection", {
            #' # Erdos-Renyi graph
            model1 <- function(n, p){
                return(igraph::sample_gnp(n, p))
            }
            # Watts-Strogatz small-world graph
            model2 <- function(n, pr, K=8){
                return(igraph::sample_smallworld(1, n, K, pr))
            }

            G <- model2(n=30, p=0.5)
            parameters <- list(seq(0.1, 0.9, 0.01), seq(0.1, 0.9, 0.01))
            result2 <- graph.model.selection(G, list(model1, model2), parameters)

            choice <- as.numeric(which.min(result2$estimate[, 2]))
            test_that("graph.model.selection is working for undirected graphs", {
                          expect_equal(choice, 2)
            })
})
