#' Distance Matrix on a List of Graphs
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
#' @param dist string indicating if you want to use the 'JS' (default), 'L1' or 'L2'
#' distances. 'JS' means Jensen-Shannon divergence.
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
#' g <- list()
#' for(i in 1:5){
#'   g[[i]] <- igraph::sample_gnp(n=50, p=0.1)
#' }
#' for(i in 6:10){
#'   g[[i]] <- igraph::sample_gnp(n=50, p=0.5)
#' }
#' for(i in 11:15){
#'  g[[i]] <- igraph::sample_gnp(n=50, p=0.9)
#' }
#' graph.dist(g, dist = 'JS')
#'
#'
#' @import stats
#' @import methods
#' @export
graph.dist <- function(Graphs, dist = "JS", symmetric = TRUE, ...) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }

    Graphs <- set.list.spectral.density(Graphs, ...)

    nGraphs <- length(Graphs)
    D <- matrix(0, nGraphs, nGraphs)

    for (i in 1:nGraphs) {
        # If the distance is symmetric, don`t bother computing both dist(x, y) and dist(y, x)
        s <- if(symmetric) { i } else { 1 }

        for (j in s:nGraphs) {
            D[i, j] <- distance(Graphs[[i]]$density, Graphs[[j]]$density, dist = dist)
        }
    }

    if(symmetric){
      D[lower.tri(D)] <- t(D)[lower.tri(D)]
    }

    return(D)
}



# distance between two density functions
distance <- function(f1, f2, dist = "KL") {
    if (dist == "KL") {
        return(KL(f1, f2))
    }
    if (dist == "L1") {
        return(L1(f1, f2))
    }
    if (dist == "L2") {
        return(L2(f1, f2))
    }
    if (dist == "JS") {
        return(JS(f1, f2))
    }
    # error if the distance parameter does not exists
    stop(paste0(dist, " distance measure is not valid. Use: KL, JS ,L1, or L2"))
}

# Return the L1 norm between two densities. This works with both real and complex distributions.
L1 <- function(f1, f2) {
    y <- abs(f1$y - f2$y)
    return(trapezoidSum(f1$x, y))
}


# Returns the L2 norm between two densities. This works with both real and complex distributions.
L2 <- function(f1, f2) {
    y <- (f1$y - f2$y)^2
    return(trapezoidSum(f1$x, y))
}

# Returns the Kullback-Leibler divergence between two densities
KL.real <- function(f1, f2) {
    y <- f1$y
    diff_zero <- (y != 0)
    eq_zero <- (f2$y == 0)

    if (sum(diff_zero & eq_zero) != 0) {
        return(Inf)
    }

    y[diff_zero] <- y[diff_zero] * (log(y[diff_zero]) - log(f2$y[diff_zero]))

    # y <- f1$y n <- length(y) for (i in 1:n) { if (y[i] != 0 && f2$y[i] == 0){ return (Inf) } if (y[i] != 0) y[i] <- y[i]*log(y[i]/f2$y[i]) }
    return(trapezoidSum.real(f1$x, y))
}

KL.complex <- function(f1, f2){
  y <- f1$y
  y <- ifelse(y != 0, y*log(y/f2$y), y)
  return(trapezoidSum.complex(f1$x, y))
}

KL <- function(f1, f2){
  if(is.matrix(f1$y)){
    KL.complex(f1, f2)
  } else {
    KL.real(f1, f2)
  }
}

# Returns the Jensen-Shannon divergence between two densities. This works with both real and complex eigenvalues (I THINK!)
JS <- function(f1, f2) {
    fm <- f1
    fm$y <- (f1$y + f2$y)/2
    return((KL(f1, fm) + KL(f2, fm))/2)
}

