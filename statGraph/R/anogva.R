#' Analysis Of Graph Variability (ANOGVA)
#'
#' \code{anogva} statistically tests whether two or more sets of graphs are generated
#' by the same random graph model. It is a generalization of the \code{takahashi.test}
#' function.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param labels an array of integers indicating the labels of each graph.
#'
#' @param maxBoot integer indicating the number of bootstrap resamplings (default \code{1000}).
#'
#' @param dist string indicating if you want to use the 'KL' (default), 'JS' , 'L1' or 'L2'
#' distances. 'KL' means Kullback-Leibler divergence. 'JS' means Jensen-Shannon divergence.
#'
#' @param ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#' @return A list with class 'htest' containing the following components:
#' \item{\code{statistic:}}{ the statistic of the test.}
#' \item{\code{p.value:}}{ the p-value of the test.}
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{data.name:}}{a string with the data's name(s).}
#'
#' @keywords analysis_of_graph_variability
#'
#' @references
#'
#' Fujita, A., Vidal, M. C. and Takahashi, D. Y. (2017) A Statistical Method to
#' Distinguish Functional Brain Networks. _Front. Neurosci._, *11*, 66.
#' doi:10.3389/fnins.2017.00066.
#'
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
#' @examples
#'
#' set.seed(1)
#' g1 <- g2 <- g3 <- list()
#' for (i in 1:20) {
#'   g1[[i]] <- igraph::sample_gnp(50, 0.50)
#'   g2[[i]] <- igraph::sample_gnp(50, 0.50)
#'   g3[[i]] <- igraph::sample_gnp(50, 0.52)
#' }
#' G <- c(g1, g2, g3)
#' label <- c(rep(1,20),rep(2,20),rep(3,20))
#' result <- anogva(G, label, maxBoot=50)
#' result
#'
#' @export
anogva <- function(Graphs, labels, maxBoot = 1000, dist = "KL", ...) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }
    data.name <- deparse(substitute(Graphs))
    # compute the spectral densities of the graphs and store in the density attribute
    Graphs <- set.list.spectral.density(Graphs, ...)
    # list of means for each diff label
    group_mean_den <- list()
    # recover unique labels
    ulabels = unique(labels)
    for (label in ulabels) {
        group_mean_den[[label]] <- get.mean.spectral.density(Graphs[labels == label])
    }
    # compute mean of means
    all_mean_den <- Graphs[[1]]$density
    all_mean_den$y <- Reduce(f = "+", Map(f = function(d) {
        d$y
    }, group_mean_den))/length(ulabels)

    distOrig <- 0
    for (label in ulabels) {
        distOrig <- distOrig + distance(group_mean_den[[label]], all_mean_den, dist = dist)
    }
    distOrig <- distOrig/length(ulabels)

    ## Permutation test
    distBoot <- array(0, maxBoot)

    for (boot in 1:maxBoot) {
        sample_labels <- sample(labels, length(labels), replace = FALSE)

        sample_group_mean_den <- list()
        # recover unique labels
        for (label in ulabels) {
            sample_group_mean_den[[label]] <- get.mean.spectral.density(Graphs[sample_labels == label])
        }

        for (label in ulabels) {
            distBoot[boot] <- distBoot[boot] + distance(sample_group_mean_den[[label]], all_mean_den, dist = dist)
        }
    }

    pvalue <- length(which(distBoot >= distOrig))/(maxBoot + 1)
    ###
    statistic <- distOrig
    names(statistic) <- "statistic"
    method_info <- paste0("Analysis of Graph Variability with\n       simulated p-value (based on ", maxBoot, " replicates)")
    rval <- list(statistic = statistic, p.value = pvalue, method = method_info, data.name = data.name)
    class(rval) <- "htest"
    return(rval)

}
