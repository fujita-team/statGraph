#' Graph spectral entropy
#'
#' \code{graph.entropy} returns the spectral entropy of an undirected graph.
#'
#' @param Graph the undirected graph (igraph object).
#' If \code{Graph} has the  attribute \code{eigenvalues} containing
#' the eigenvalues of \code{Graph}, such values will be used to
#' compute its spectral density.
#'
#' @param ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#' @return a real number corresponding to the graph spectral entropy.
#'
#' @keywords spectral_entropy
#'
#' @references
#' Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
#' Discriminating Different Classes of Biological Networks by Analyzing the
#' Graph Spectra  Distribution. _PLoS ONE_, *7*, e49949.
#' doi:10.1371/journal.pone.0049949.
#'
#' @examples
#' G <- igraph::sample_gnp(n=100, p=0.5)
#' entropy <- graph.entropy(Graph = G)
#' entropy
#'
#' @export
graph.entropy <- function(Graph,...) {
  if(!valid.input(Graph)) stop("The input should be an igraph object!")
  # compute the entropy
  f <- graph.spectral.density(Graph, ...)

  if (sum(is.na(f)) > 0)
    return(NA)
  y <- f$y
  valid_idx <- which(y != 0)
  y[valid_idx] <- y[valid_idx]*log(y[valid_idx])

  return(-trapezoidSum(f$x, y))
}
