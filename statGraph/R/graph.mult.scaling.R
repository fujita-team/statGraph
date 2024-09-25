#' Multidimensional Scaling of Graphs
#'
#' \code{graph.mult.scaling} performs multidimensional scaling of graphs. It
#' takes the Jensen-Shannon divergence between graphs (JS) and uses the
#' \code{cmdscale} function from the \code{stats} package to obtain a set of points such
#' that the distances between the points are similar to JS.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param plot logical. If \code{TRUE} (default) the points chosen to represent the
#' Jensen-Shannon divergence between graphs are plotted.
#'
#' @param type what type of plot should be drawn. The default value is \code{'n'},
#' which indicates that the points will not be plotted (i. e. only the labels
#' of the graphs will be plotted).
#'
#' @param dist string indicating if you want to use the 'JS' (default), 'L1' or 'L2'
#' distances. 'JS' means Jensen-Shannon divergence.
#'
#' @param main title of the plot (default value is '').
#'
#' @param ... additional parameters for \code{\link{graph.spectral.density}}.
#'
#' @return A list with class 'statGraph' containing the following components:
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{info:}}{ a string showing details about the method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' \item{\code{values:}}{ a matrix in which each column corresponds to a coordinate and each
#' row corresponds to a graph (point). Then, each row gives the coordinates of
#' the points chosen to represent the Jensen-Shannon divergence (by default), L1, or
#' L2 distance between graphs.}
#'
#' @keywords multidimensional_scaling
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
#' graph.mult.scaling(G)
#'
#' @importFrom graphics text
#' @export
graph.mult.scaling <- function(Graphs, plot = TRUE, type = "n", dist = "JS", main = "", ...) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }

    data.name <- deparse(substitute(Graphs))

    Graphs <- set.list.spectral.density(Graphs, ...)
    nGraphs <- length(Graphs)

    d <- graph.dist(Graphs = Graphs, dist = dist, ...)

    if (is.null(names(Graphs)))
        names <- as.character(1:nGraphs) else names <- names(Graphs)
    colnames(d) <- rownames(d) <- names
    fit <- cmdscale(as.dist(d), k = 2)

    x <- fit[, 1]
    y <- fit[, 2]
    names(x) <- names
    names(y) <- names
    if (plot) {
        plot(x, y, xlab = "x", ylab = "y", main = main, type = type)
        text(x, y, labels = names, cex = 1)
    }
    method <- "Multidimensional scaling of graphs"
    info <- ""
    output <- list(method = method, info = info, values = fit)
    class(output) <- "statGraph"
    return(output)
}
