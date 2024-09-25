#' K-means for Graphs
#'
#' \code{graph.kmeans} clusters graphs following a k-means algorithm based on the
#' Jensen-Shannon divergence between the spectral densities of the graphs.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param k an integer specifying the number of clusters.
#'
#' @param nstart the number of trials of k-means clusterizations. The algorithm
#' returns the clusterization with the best silhouette.
#'
#' @param dist string indicating if you want to use the 'JS' (default), 'L1' or 'L2'
#' distances. 'JS' means Jensen-Shannon divergence.
#'
#' @param  ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#' @return A list with class 'statGraph' containing the following components:
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{info:}}{ a string showing details about the method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' \item{\code{cluster:}}{ a vector of the same length of \code{Graphs} containing the clusterization
#' labels.}
#' \item{\code{centers:}}{a list containing the centroids of each cluster.}
#'
#' @keywords k-means
#'
#' @references
#' MacQueen, James. 'Some methods for classification and analysis of
#' multivariate observations.' Proceedings of the fifth Berkeley symposium on
#' mathematical statistics and probability. Vol. 1. No. 14. 1967.
#'
#' Lloyd, Stuart. 'Least squares quantization in PCM.' IEEE transactions on
#' information theory 28.2 (1982): 129-137.
#'
#' @examples
#' set.seed(1)
#' g <- list()
#' for(i in 1:5){
#'   g[[i]] <- igraph::sample_gnp(30, p=0.2)
#' }
#' for(i in 6:10){
#'   g[[i]] <- igraph::sample_gnp(30, p=0.5)
#' }
#' res <- graph.kmeans(g, k=2, nstart=2)
#' res
#'
#' @export
graph.kmeans <- function(Graphs, k, nstart = 2, dist = "JS", ...) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }

    data.name <- deparse(substitute(Graphs))

    nGraphs <- length(Graphs)
    Graphs <- set.list.spectral.density(Graphs, ...)

    sil <- -1

    if (k > nstart) {
        nstart <- k
    }
    label.final <- NULL
    centroid.final <- NULL
    for (ns in 1:nstart) {
        ## random initialization of the clusters
        label <- sample(rep_len(1:k, nGraphs))

        converged <- FALSE
        while (converged == FALSE) {
            centroid <- list()

            for (j in 1:k) {
                centroid[[j]] <- get.mean.spectral.density(Graphs[label == j])
            }

            distance_mat <- matrix(0, nGraphs, k)
            for (j in 1:k) {
                for (i in 1:nGraphs) {
                  distance_mat[i, j] <- distance(Graphs[[i]]$density, centroid[[j]], dist = dist)
                }
            }

            label.new <- array(0, nGraphs)
            for (i in 1:nGraphs) {
                label.new[i] <- which(distance_mat[i, ] == min(distance_mat[i, ]))[1]
            }
            i <- 1
            while (i <= k) {
                if (length(which(label.new == i)) != 0) {
                  i <- i + 1
                } else {
                  ## there is an empty cluster
                  size.cluster <- array(0, k)
                  for (j in 1:k) {
                    size.cluster[j] <- length(which(label.new == j))
                  }
                  largest.cluster <- which(size.cluster == max(size.cluster))
                  item <- which(distance_mat[, largest.cluster] == max(distance_mat[which(label.new == largest.cluster), largest.cluster]))
                  label.new[item] <- i
                  i <- 1
                }
            }

            if (sum(label == label.new) == nGraphs) {
                converged <- TRUE
                sil.new <- mean(cluster::silhouette(label, graph.dist(Graphs, dist = dist))[, 3])

                if (sil.new > sil) {
                  sil <- sil.new
                  label.final <- label
                  centroid.final <- centroid
                }
            }
            label <- label.new
        }
    }
    ###
    method_info <- "K-means for Graphs"
    info <- "Clustering the graphs following a k-means algorithm"
    output <- list(method = method_info, info = info, data.name = data.name, cluster = label.final, centers = centroid.final)
    class(output) <- "statGraph"
    return(output)
}
