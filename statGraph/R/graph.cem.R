#' Graph Clustering Expectation-Maximization (gCEM)
#'
#' \code{graph.cem} clusters graphs following an expectation-maximization algorithm based
#' on the Kullback-Leibler divergence between the spectral densities of the
#' graph and of the random graph model.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param  model a string that indicates one of the following random graph
#' models: 'ER' (Erdos-Renyi random graph), 'GRG' (geometric random graph), 'KR'
#' (k regular graph), 'WS' (Watts-Strogatz model), and 'BA' (Barabási-Albert
#' model).
#'
#' @param k an integer specifying the number of clusters.
#'
#' @param max_iter the maximum number of expectation-maximization steps to execute.
#'
#' @param ... Other relevant parameters for \code{\link{graph.param.estimator}}.
#'
#' @return A list with class 'statGraph' containing the following components:
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{info:}}{ a string showing details about the method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' \item{\code{cluster:}}{ a vector of the same length of \code{Graphs} containing the clusterization
#' labels.}
#' \item{\code{parameters:}}{ a vector containing the estimated parameters for the groups.
#' It has the length equals to \code{k}.}
#'
#' @keywords graph.cem
#'
#' @references
#' Celeux, Gilles, and Gerard Govaert. 'Gaussian parsimonious clustering
#' models.' Pattern recognition 28.5 (1995): 781-793.
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' _Journal of the Royal Statistical Society series B_, 53, 683-690.
#' http://www.jstor.org/stable/2345597.
#'
#' @examples
#' \donttest{
#'  set.seed(1)
#'  g <- list()
#'  for(i in 1:2){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=0.5)
#'  }
#'  for(i in 3:4){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=1)
#'  }
#'  res <- graph.cem(g, model='ER', k=2, max_iter=1,eps=0.1)
#'  res
#'  }
#' @export
graph.cem <- function(Graphs, model, k, max_iter = 10, ...) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }
    data.name <- deparse(substitute(Graphs))
    ## Pre-processing of the graph spectra
    Graphs <- set.list.spectral.density(Graphs, ...)

    nGraphs <- length(Graphs)
    tau <- matrix(0, nrow = k, ncol = nGraphs)
    kl <- matrix(0, nrow = k, ncol = nGraphs)

    prevlik <- 0
    lik <- 1
    count <- 0
    prevlabels <- array(0, nGraphs)
    labels <- array(0, nGraphs)
    g_GIC <- array(0, nGraphs)
    p <- array(0, k)

    p_graph <- array(0, nGraphs)
    ## Parameter estimation
    ret <- Map(f = function(G) {
        graph.param.estimator(G, model = model, ...)
    }, Graphs)
    #
    for (i in 1:nGraphs) {
        p_graph[i] <- ret[[i]]$param
        g_GIC[i] <- ret[[i]]$dist
    }

    # Initialize cluster parameters
    p_uniq <- unique(p_graph)
    for (i in 1:k) {
        p[i] <- quantile(p_uniq, i/(k + 1))
        # the KR parameter needs to be even
        if (model == "KR") {
            p[i] <- round(p[i])
        }
    }

    converged <- 0
    count <- 0
    while (!converged) {
        kl <- matrix(0, nrow = k, ncol = nGraphs)
        for (i in 1:k) {
            for (j in 1:nGraphs) {
                kl[i, j] = GIC(Graph = Graphs[[j]], model = model, p = p[i], ...)$value
            }
        }

        kl[which(kl == Inf)] <- max(kl[which(kl < Inf)])
        kl[which(kl == 0)] <- 1e-09

        for (i in 1:nGraphs) {
            tau[, i] <- (1/kl[, i])/sum(1/kl[, i])
        }


        for (i in 1:nGraphs) {
            labels[i] <- which(tau[, i] == max(tau[, i]))[1]
        }
        # Check if there is an empty group
        for (i in 1:k) {
            if (length(which(labels == i)) == 0) {
                labels[which(tau[i, ] == max(tau[i, ]))] <- i
            }
        }

        # Estimates the value of p for the models to maximize O tae
        for (i in 1:k) {
            p[i] <- sum(p_graph[which(labels == i)])/length(which(labels == i))
            if (model == "KR") {
                p[i] <- round(p[i])
            }
        }

        prevlik <- lik
        lik <- sum(tau * kl)
        count <- count + 1
        if (count > max_iter) {
            converged <- TRUE
        }

        if ((prevlik != 0 && prevlik/lik > 0.99 && prevlik/lik < 1.01)) {
            converged <- TRUE
        }
        prevlabels <- labels
    }
    ###
    method_info <- "Graph Clustering Expectation-Maximization "
    info <- ret[[1]]$info  # info of the parameter estimator
    result <- list(cluster = labels, parameters = p)
    output <- list(method = method_info, info = info, data.name = data.name, cluster = result$cluster, parameters = result$parameters)
    class(output) <- "statGraph"
    return(output)
}




