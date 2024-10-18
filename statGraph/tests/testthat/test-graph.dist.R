RunTest( {
  A <- list()
  for(i in 1:10){
    p <- sample(seq(0.1, 0.9, 0.1), 1)
    A[[i]] <- igraph::sample_gnp(50, p)
  }

  # JS is symmetric, so d1 should be equal to d2
  d1 <- statGraph::graph.dist(A, dist="JS", symmetric = FALSE)
  d2 <- statGraph::graph.dist(A, dist="JS", symmetric = TRUE)

  test_that("graph.dist is working well", {
    expect_true(identical(d1, d2))
  })
})

RunTest( {
  A <- list()
  for(i in 1:10){
    p <- sample(seq(0.1, 0.9, 0.1), 1)
    A[[i]] <- igraph::sample_gnp(50, p)
  }

  # L1 is symmetric, so d1 should be equal to d2
  d1 <- statGraph::graph.dist(A, dist="L1", symmetric = FALSE)
  d2 <- statGraph::graph.dist(A, dist="L1", symmetric = TRUE)

  test_that("graph.dist is working well", {
    expect_true(identical(d1, d2))
  })
})

# Directed tests

RunTest( {
  A <- list()
  for(i in 1:10){
    p <- sample(seq(0.1, 0.9, 0.1), 1)
    A[[i]] <- igraph::sample_gnp(50, p, directed = TRUE)
  }

  # L2 is symmetric, so d1 should be equal to d2
  d1 <- statGraph::graph.dist(A, dist="L2", symmetric = FALSE)
  d2 <- statGraph::graph.dist(A, dist="L2", symmetric = TRUE)

  test_that("graph.dist is working well", {
    expect_true(identical(d1, d2))
  })
})

RunTest( {
  A <- list()
  for(i in 1:10){
    p <- sample(seq(0.1, 0.9, 0.1), 1)
    A[[i]] <- igraph::sample_gnp(50, p, directed = TRUE)
  }

  # JS is symmetric, so d1 should be equal to d2
  d1 <- statGraph::graph.dist(A, dist="JS", symmetric = FALSE)
  d2 <- statGraph::graph.dist(A, dist="JS", symmetric = TRUE)

  test_that("graph.dist is working well", {
    expect_true(identical(d1, d2))
  })
})


