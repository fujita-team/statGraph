# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

RunTest( {
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

RunTest( {
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


RunTest( {
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

RunTest( {
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

RunTest( {

  # Extended models
  PAE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_pa(n=N, power=1 * P + 0.001, m=N/10)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  WSE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_smallworld(dim=1, size=N, nei=N/10, p=P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  ERE <- function(n, p){
    f <- function(N, P) igraph::sample_gnp(N, P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  # A graph
  G <- ERE(100, 0.3)

  result <- statGraph::graph.model.selection(G, models = list(PAE, WSE, ERE), parameters = list(
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9)),
    , model_names = c("PAE", "WSE", "ERE"), ngraphs=5, npoints=512)

  test_that("graph.model.selection is working for directed graphs", {
    expect_equal(result$model, "ERE")
  })
})

RunTest( {

  # Extended models
  PAE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_pa(n=N, power=1 * P + 0.001, m=N/10)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  WSE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_smallworld(dim=1, size=N, nei=N/10, p=P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  ERE <- function(n, p){
    f <- function(N, P) igraph::sample_gnp(N, P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  # A graph
  G <- PAE(100, 0.3)

  result <- statGraph::graph.model.selection(G, models = list(PAE, WSE, ERE), parameters = list(
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9)),
    , model_names = c("PAE", "WSE", "ERE"), ngraphs=5, npoints=512)

  test_that("graph.model.selection is working for directed graphs", {
    expect_equal(result$model, "PAE")
  })
})

RunTest( {

  # Extended models
  PAE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_pa(n=N, power=1 * P + 0.001, m=N/10)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  WSE <- function(n, p)
  {
    f <- function(N, P) igraph::sample_smallworld(dim=1, size=N, nei=N/10, p=P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  ERE <- function(n, p){
    f <- function(N, P) igraph::sample_gnp(N, P)
    statGraph:::graph.extended.model(f, n, p, p)
  }

  # A graph
  G <- WSE(100, 0.3)

  result <- statGraph::graph.model.selection(G, models = list(PAE, WSE, ERE), parameters = list(
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9),
    list(lo=0.1, hi=0.9)),
    , model_names = c("PAE", "WSE", "ERE"), ngraphs=5, npoints=512)

  test_that("graph.model.selection is working for directed graphs", {
    expect_equal(result$model, "WSE")
  })
})
