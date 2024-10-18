# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf


RunTest({
    Graphs <- list()
    labels <- c()

    for(i in 1:10){
      G1 <- igraph::sample_gnp(100, 0.3, directed = FALSE)
      G2 <- igraph::sample_gnp(100, 0.5, directed = FALSE)
      Graphs <- c(Graphs, list(G1, G2))
      labels <- c(labels, c(1, 2))
    }

    result <- statGraph::graph.kmeans(Graphs, 2)
    ari <- mclust::adjustedRandIndex(result$cluster, labels)
    test_that("graph.kmeans works with undirected graphs.", {
      expect_equal(ari, 1)
    })
})

RunTest({
  Graphs <- list()
  labels <- c()

  for(i in 1:10){
    G1 <- igraph::sample_gnp(100, 0.3, directed = TRUE)
    G2 <- igraph::sample_gnp(100, 0.5, directed = TRUE)
    Graphs <- c(Graphs, list(G1, G2))
    labels <- c(labels, c(1, 2))
  }

  result <- statGraph::graph.kmeans(Graphs, 2)
  ari <- mclust::adjustedRandIndex(result$cluster, labels)
  test_that("graph.kmeans works with directed graphs.", {
    expect_equal(ari, 1)
  })
})
