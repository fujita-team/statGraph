#' Test for Association / Correlation Between Paired Samples of Graphs
#'
#' \code{graph.cor.test} tests for association between paired samples of graphs,
#' using Spearman's rho correlation coefficient.
#'
#' @param Graphs1 a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param Graphs2 a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @return A list with class 'htest' containing the following components:
#' \item{\code{statistic:}}{ the value of the test statistic.}
#' \item{\code{p.value:}}{ the p-value of the test.}
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' \item{\code{estimates:}}{ the estimated measure of association 'rho'.}
#'
#' @keywords correlation_coefficient
#'
#' @references
#' Fujita, A., Takahashi, D. Y., Balardin, J. B., Vidal, M. C. and Sato, J. R.
#' (2017) Correlation between graphs with an application to brain network
#' analysis. _Computational Statistics & Data Analysis_ *109*, 76-92.
#'
#' @examples
#' library(mvtnorm)
#'
#' set.seed(1)
#' G1 <- G2 <- list()
#'
#' p <- mvtnorm::rmvnorm(50, mean=c(0,0), sigma=matrix(c(1, 0.5, 0.5, 1), 2, 2))
#'
#' ma <- max(p)
#' mi <- min(p)
#' p[,1] <- (p[,1] - mi)/(ma - mi)
#' p[,2] <- (p[,2] - mi)/(ma - mi)
#'
#' for (i in 1:50) {
#'   G1[[i]] <- igraph::sample_gnp(50, p[i,1])
#'   G2[[i]] <- igraph::sample_gnp(50, p[i,2])
#' }
#' graph.cor.test(G1, G2)
#'
#' @import stats
#' @import methods
#' @import MASS
#' @import mvtnorm
#' @export
#'
graph.cor.test <- function(Graphs1, Graphs2) {
    if (!valid.input(Graphs1, level = 1) || !valid.input(Graphs2, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }

    data.name <- paste(deparse(substitute(Graphs1)), "and", deparse(substitute(Graphs2)))

    G1.radius <- unlist(Map(f = function(G) {
        get.largest.eigenvalue(G)
    }, Graphs1))
    G2.radius <- unlist(Map(f = function(G) {
        get.largest.eigenvalue(G)
    }, Graphs2))

    res <- cor.test(G1.radius, G2.radius, method = "spearman")
    ###
    statistic <- res$statistic
    names(statistic) <- "statistic"
    estimate <- res$estimate
    names(estimate) <- "rho"
    method <- "Association between paired samples of graphs, using Spearman's rho correlation coefficient"
    rval <- list(statistic = statistic, p.value = res$p.value, method = method, data.name = data.name, estimate = estimate)
    class(rval) <- "htest"
    return(rval)
}
