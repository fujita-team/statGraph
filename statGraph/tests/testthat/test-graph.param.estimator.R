# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

RunTest("graph.param.estimator", {
    tolerance <- 0.01
    model <- function(n, p) {
        return(igraph::sample_gnp(n, p))
    }
    G <- model(100, 0.3)
    result <- graph.param.estimator(G, model,  seq(0.2, 0.8, 0.1))
    test_that("Documentation example is working",
              {
                expect_equal(result$param, 0.3, tolerance = tolerance)
              }
    )
})

RunTest("graph.param.estimator", {
    tolerance <- 0.01
    model <- function(n, p) {
        return(igraph::sample_gnp(n, p))
    }
    G <- model(100, 0.3)
    result <- graph.param.estimator(G, model,  seq(0.2, 0.8, 0.1), search="ternary")
    test_that("Ternary search is working",
              {
                expect_equal(result$param, 0.3, tolerance = tolerance)
              }
    )
})

RunTest("graph.param.estimator", {
    tolerance <- 0.01
    model <- function(n, p) {
        return(igraph::sample_gnp(n, p))
    }
    G <- model(100, 0.3)
    result <- graph.param.estimator(G, model,  seq(0.2, 0.8, 0.1), dist = "L1")
    test_that("Grid search with L1",
              {
                expect_equal(result$param, 0.3, tolerance = tolerance)
              }
    )
})

RunTest("graph.param.estimator", {
    tolerance <- 0.01
    model <- function(n, p) {
        return(igraph::sample_gnp(n, p))
    }
    G <- model(100, 0.3)
    result <- graph.param.estimator(G, model,  seq(0.2, 0.8, 0.1), search="ternary", dist = "L2")
    test_that("Ternary search with L1",
              {
                expect_equal(result$param, 0.3, tolerance = tolerance)
              }
    )
})
