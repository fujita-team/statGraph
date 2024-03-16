#' Distance matrix on a list of graphs.
#'
#' Given a list of graphs, \code{graph.dist} builds a distance matrix
#' according to the Jensen-Shannon divergence, L2 norm, or L1 norm
#' between the spectral density of the graphs graphs.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param dist string indicating if you want to use the "JS" (default), "L1" or "L2"
#' distances. "JS" means Jensen-Shannon divergence.
#'
#' @param ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#'
#' @return a distance matrix
#'
#' @keywords distance_matrix
#'
#' @examples
#' set.seed(1)
#' G <- list()
#' for (i in 1:5) {
#'   G[[i]] <- igraph::sample_gnp(50, 0.5)
#' }
#' graph.dist(G, 3)
#'
#' @import stats
#' @import methods
#' @export
graph.dist <- function(Graphs,dist = "JS",...) {
  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")

  Graphs <- set.list.spectral.density(Graphs,...)

  nGraphs <- length(Graphs)
  D <- matrix(0, nGraphs, nGraphs)
  for (i in 1:(nGraphs-1)) {
    for (j in (i+1):nGraphs) {
      D[i,j] <- D[j,i] <- distance(Graphs[[i]]$density,Graphs[[j]]$density,dist = dist)
    }
  }
  D = as.dist(D)
  return(D)
}
