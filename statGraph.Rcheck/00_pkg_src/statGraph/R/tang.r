#' Tang hypothesis testing for random graphs.
#'
#' Given two independent finite-dimensional random dot product graphs,
#' \code{tang.test} tests if they have generating latent positions that are drawn
#' from the same distribution.
#'
#' @param Graph1 the first undirected graph (igraph object).
#' If \code{Graph} has the  attribute \code{eigenvalues} containing
#' the eigenvalues of \code{Graph}, such values will be used to
#' compute its spectral density.
#'
#' @param Graph2 the second undirected graph (igraph object).
#' If \code{Graph} has the  attribute \code{eigenvalues} containing
#' the eigenvalues of \code{Graph}, such values will be used to
#' compute its spectral density.
#'
#' @param dim dimension of the adjacency spectral embedding.
#'
#' @param sigma a real value indicating the kernel bandwidth. If \code{NULL} (default)
#' the bandwidth is calculated by the method.
#'
#' @param alpha the significance level for the test (default is \code{0.05}).
#'
#' @param bootstrap_sample integer indicating the number of bootstrap resamples
#' (default is \code{200}).
#'
#' @param printResult logical indicating if the test must print the result
#' (default is \code{FALSE}).
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{X1:}}{ the embedding of \code{Graph1}.}
#' \item{\code{X2:}}{ the embedding of \code{Graph2}.}
#' \item{\code{test_stats:}}{ the value of the test.}
#' \item{\code{p_value:}}{ the p-value of the test.}
#' \item{\code{bootstrap_samples:}}{ the test distribution on the bootstrap resamples.}
#' }
#'
#' @references
#' Tang, Minh, et al. "A nonparametric two-sample hypothesis testing problem for
#' random graphs." Bernoulli 23.3 (2017): 1599-1630.
#'
#' Tang, Minh, et al. "A semiparametric two-sample hypothesis testing problem
#' for random graphs." Journal of Computational and Graphical Statistics 26.2
#' (2017): 344-354.
#'
#' @examples
#' set.seed(42)
#'
#' ## test under H0
#' lpvs <- matrix(rnorm(200), 20, 10)
#' lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
#' Graph1 <- igraph::sample_dot_product(lpvs)
#' Graph2 <- igraph::sample_dot_product(lpvs)
#' D <- graph.tang.test(Graph1,Graph2, 5, printResult = TRUE)
#'
#' ## test under H1
#' lpvs2 <- matrix(pnorm(200), 20, 10)
#' lpvs2 <- apply(lpvs2, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
#' g2 <- suppressWarnings(igraph::sample_dot_product(lpvs2))
#' D <- graph.tang.test(Graph1,Graph2, 5, printResult = TRUE)
#'
#'
#' @export
graph.tang.test <- function(Graph1, Graph2, dim, sigma = NULL, alpha = 0.05, bootstrap_sample=200,
                            printResult = FALSE){

  tang.validateInput(Graph1, Graph2, dim, alpha, bootstrap_sample,printResult)
  Xhat1 = tang.embed.graph(Graph1, dim)
  Xhat2 = tang.embed.graph(Graph2, dim)
  if(is.null(sigma)){
    sigma = tang.get.sigma(Xhat1, Xhat2)
  }
  test_stat = tang.test.stat(Xhat1, Xhat2, sigma)
  test_distribution = tang.sampling.distribution(Graph1, dim, bootstrap_sample)
  #test_distribution2 = sampling.distribution(G2, dim, bootstrap_sample)
  p_val = tang.pvalue(test_stat, test_distribution)
  num = round(bootstrap_sample * alpha)
  reject = FALSE
  if (p_val <= alpha) {
    reject = TRUE
  }
  if (printResult) {
    if (reject) {
      msg <- paste("Reject the null hypothesis that two graphs are identically",
                   "distributed.")
      print(msg)
    } else {
      msg <- paste("Fail to reject the null hypothesis that two graphs are",
                   "identically distributed.")
      print(msg)
    }
  }
  out = list(X1 = Xhat1, X2 = Xhat2, test_stats = test_stat,
             p_value = p_val, bootstrap_samples = test_distribution)
  return(out)
}


## Auxiliary for Tang method.
tang.test.stat <- function(X, Y, sigma) {
  n <- nrow(X)
  m <- nrow(Y)
  tmpXX <- sum(exp(-(as.matrix(stats::dist(X))^2)/(2*sigma^2)))
  tmpYY <- sum(exp(-(as.matrix(stats::dist(Y))^2)/(2*sigma^2)))
  tmpXY <- sum(exp(-(tang.rect.dist(X,Y))/(2*sigma^2)))
  tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
  return((m+n)*tmp)
}

## Auxiliary for Tang method.
tang.embed.graph <- function(Graph, dim) {
  # Call ase to get latent position
  default_options = igraph::arpack_defaults()
  default_options$maxiter = .Machine$integer.max
  lpv = igraph::embed_adjacency_matrix(Graph,dim, options = default_options)$X

  # Fix signs of eigenvectors issue
  for (i in 1:dim) {
    if (sign(lpv[1, i]) != 1) {
      lpv[, i] = -lpv[, i]
    }
  }
  return(lpv)
}

## Auxiliary for Tang method.
tang.rect.dist <- function(X,Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))

  D <- tmp2 - 2*tmp1 + tmp3
  return(D)
}

## Auxiliary for Tang method.
tang.get.sigma <- function(X1, X2) {
  v1 = as.vector(stats::dist(X1))
  v2 = as.vector(stats::dist(X2))
  v = base::append(v1, v2)
  sigma = stats::median(v)
  return(sigma)
}

## Auxiliary for Tang method.
tang.sampling.distribution <- function(G1, dim, bootstrap_sample_size) {
  Xhat1 = tang.embed.graph(G1,dim)
  P = t(Xhat1)
  test_distribution = c()
  i = 1
  while (i <= bootstrap_sample_size) {
    tryCatch({
      G_a = suppressWarnings(igraph::sample_dot_product(P))
      G_b = suppressWarnings(igraph::sample_dot_product(P))
      Xhat_a = suppressWarnings(tang.embed.graph(G_a, dim))
      Xhat_b = suppressWarnings(tang.embed.graph(G_b, dim))
      sigma = tang.get.sigma(Xhat_a, Xhat_b)
      ts = tang.test.stat(Xhat_a, Xhat_b, sigma)
      test_distribution[i] = ts
      i = i + 1
    }, error=function(e) {stop(print(e))})
  }
  test_distribution
}

## Auxiliary for Tang method.
tang.pvalue <- function(ts, test_distribution) {
  area = sum(test_distribution >= ts) / length(test_distribution)
  return(area)
}

## Auxiliary for Tang method.
tang.validateInput <- function(G1, G2, dim, alpha, bootstrap_sample, printResult) {

  if (!methods::is(G1,'igraph')) { stop("Input object 'Graph1' is not an igraph object.") }
  if (!methods::is(G2,'igraph')) { stop("Input object 'Graph2' is not an igraph object.") }
  if (!is.null(dim)) {
    if (!methods::is(dim,"numeric")) { stop("Input 'dim' is not a number.") }
    if (dim%%1 != 0) { stop("Input 'dim' must be an integer.") }
    if (dim < 1) { stop("Number of dimensions 'dim' is less than 1.") }
    if (dim >= igraph::vcount(G1) || dim >= igraph::vcount(G2)) { stop("Num. Embedded dimensions 'dim' is greater or equal than number of vertices.") }
  }

  if (!methods::is(alpha,"numeric")) {
    stop("Input object 'alpha' is not a numeric value.")
  } else {
    if (alpha >= 1 || alpha <= 0) {
      stop("Significance level alpha must be strictly between 0 and 1.")
    }
  }
  if (!methods::is(bootstrap_sample,"numeric")) {
    stop("Input object 'bootstrap_sample' is not a numeric value.")
  } else {
    if (bootstrap_sample <= 20) {
      stop("The size of bootstrap sample is too small. Pick a larger value.")
    }
  }
  if (!is.logical(printResult)) { stop("Error: Input 'printResult' must be a logical.")}
}

