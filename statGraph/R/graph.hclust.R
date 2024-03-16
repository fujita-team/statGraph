#' Hierarchical cluster analysis on a list of graphs.
#'
#' Given a list of graphs, \code{graph.hclust} builds a hierarchy of clusters
#' according to the Jensen-Shannon divergence between graphs.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param k the number of clusters.
#'
#' @param clus_method the agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of '"ward.D"', '"ward.D2"', '"single"',
#' '"complete"', '"average"' (= UPGMA), '"mcquitty"' (= WPGMA), '"median"'
#' (= WPGMC) or '"centroid"' (= UPGMC).
#'
#' @param dist string indicating if you want to use the "JS" (default), "L1" or "L2"
#' distances. "JS" means Jensen-Shannon divergence.
#'
#' @param ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{hclust:}}{ an object of class \code{hclust} which describes the tree produced
#' by the clustering process.}
#' \item{\code{cluster:}}{ the clustering labels for each graph.}
#' }
#'
#'
#' @keywords clustering
#'
#' @references
#' Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
#' Discriminating Different Classes of Biological Networks by Analyzing the
#' Graph Spectra  Distribution. _PLoS ONE_, *7*, e49949.
#' doi:10.1371/journal.pone.0049949.
#'
#' Silverman, B. W. (1986) _Density Estimation_.  London: Chapman and Hall.
#'
#' Sturges, H. A. The Choice of a Class Interval. _J. Am. Statist. Assoc._,
#' *21*, 65-66.
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' _Journal of the Royal Statistical Society series B_, 53, 683-690.
#' http://www.jstor.org/stable/2345597.
#'
#' @examples
#' set.seed(1)
#' G <- list()
#' for (i in 1:5) {
#'   G[[i]] <- igraph::sample_gnp(50, 0.5)
#' }
#' for (i in 6:10) {
#'   G[[i]] <- igraph::sample_smallworld(1, 50, 8, 0.2)
#' }
#' for (i in 11:15) {
#'   G[[i]] <- igraph::sample_pa(50, power = 1, directed = FALSE)
#' }
#' graph.hclust(G, 3)
#'
#' @import stats
#' @import methods
#' @export
graph.hclust <- function(Graphs, k, clus_method="complete", dist = "JS", ...) {

  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")

  Graphs <- set.list.spectral.density(Graphs,...)

  d <- graph.dist(Graphs,dist = dist)

  tmp <- hclust(d, method=clus_method)

  res <- list()
  res$hclust <- tmp
  res$cluster <- cutree(tmp, k)

  return(res)
}
