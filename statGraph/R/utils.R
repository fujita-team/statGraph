# obtain the eigenvalues of the graph, if Graph contains eigenvalues as attribute then such values are returned
graph.eigenvalues <- function(Graph) {
    if (!is.null(Graph$density)) {
        return(NULL)
    }
    if (!is.null(Graph$eigenvalues)) {
        return(Graph$eigenvalues)
    } else {
        A <- igraph::as_adjacency_matrix(Graph)
        eigenvalues <- as.numeric(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
        eigenvalues <- eigenvalues/sqrt(nrow(A))
        rm(A)
        return(eigenvalues)
    }
}

# checks if the variable is a graph, a list of graphs, or a list of lists of graphs level = 0: checks if input is a graph level = 1: checks if the input is a list of graphs level = 2: checks if the input is a list of list of
# graphs, and so on
valid.input <- function(Graph, level = 0) {
    if (level == 0) {
        return(methods::is(Graph, "igraph"))
    } else {
        if (methods::is(Graph, "list")) {
            return(Reduce(f = "&", Map(f = function(x) {
                valid.input(x, level - 1)
            }, Graph)))
        } else {
            return(FALSE)
        }
    }
}


# Returns the density function for a sample x at n points in the interval [from, to]
gaussianDensity <- function(x, from = NULL, to = NULL, bandwidth = "Silverman", npoints = 1024) {
    # if all values of x are the same, only SIlverman works
    if ((max(x) == min(x) && (bandwidth != "Silverman"))) {
        # this case happens when all eigenvalues are equal, and the used bandwidth is sturges
        stop(paste0("All eigenvalues are equal to ", round(x[1], 3), " use 'Silverman' bandwidth instead."))
    }
    #
    if (bandwidth == "Sturges") {
        bw <- kernelBandwidth(x)
    } else if (bandwidth == "Silverman") {
        bw <- stats::bw.nrd0(x)
    } else if (bandwidth == "bcv") {
        bw <- suppressWarnings(stats::bw.bcv(x))
    } else if (bandwidth == "ucv") {
        bw <- suppressWarnings(stats::bw.ucv(x))
    } else if (bandwidth == "SJ") {
        bw <- "SJ"
    } else {
        stop("Please, choose a valid bandwidth.")
    }

    if (is.null(from) || is.null(to)) {
        f <- stats::density(x, bw = bw, n = npoints)
    } else {
        f <- stats::density(x, bw = bw, from = from, to = to, n = npoints)
    }
    f$y <- f$y + 1e-12  # we do not want the area to be zero, so we add a very small number
    area <- trapezoidSum(f$x, f$y)
    return(list(x = f$x, y = f$y/area, from = min(f$x), to = max(f$x), bw = f$bw, method = "exact"))
}


# Given a partition x[1]...x[n] and y[i] = f(x[i]), returns the trapezoid sum approximation for int_{x[1]}^{x[n]}{f(x)dx}
trapezoidSum <- function(x, y) {
    n <- length(x)
    delta <- (x[2] - x[1])
    area <- sum(y[2:(n - 1)])
    area <- (area + (y[1] + y[n])/2) * delta
    return(area)
}


# Returns the kernel bandwidth for a sample x based on Sturge's criterion
kernelBandwidth <- function(x) {
    n <- length(x)
    nbins <- ceiling(log2(n) + 1)
    return(abs(max(x) - min(x))/nbins)
}

# functions to obtain the smallest and largest eigenvalues of a graph or a list of graphs
get.smallest.eigenvalue <- function(Graphs) {
    if (methods::is(Graphs, "igraph")) {
        if (is.null(Graphs$eigenvalues)) {
            A <- igraph::as_adjacency_matrix(Graphs, type = "both")
            ev <- rARPACK::eigs_sym(A, k = 1, which = "SA")$values[1]
            rm(A)
            return(ev)
        } else {
            return(Graphs$eigenvalues[igraph::vcount(Graphs)])
        }
    } else if (methods::is(Graphs, "list")) {
        return(Reduce(f = "min", Map(f = get.smallest.eigenvalue, Graphs)))
    }
    stop("Input should be a Graph or a list of graphs.")
}

# functions to obtain the largest and largest eigenvalues of a graph or a list of graphs
get.largest.eigenvalue <- function(Graphs) {
    if (methods::is(Graphs, "igraph")) {
        if (is.null(Graphs$eigenvalues)) {
            A <- igraph::as_adjacency_matrix(Graphs, type = "both")
            ev <- rARPACK::eigs_sym(A, k = 1)$values[1]
            rm(A)
            return(ev)
        } else {
            return(Graphs$eigenvalues[1])
        }
    } else if (methods::is(Graphs, "list")) {
        return(Reduce(f = "max", Map(f = get.largest.eigenvalue, Graphs)))
    }
    stop("Input should be a Graph or a list of graphs.")
}

# transform a list of graphs to their respective adjacency matrices
graphList.to.adjList <- function(Graphs) {
    if (methods::is(Graphs, "igraph")) {
        return(igraph::as_adjacency_matrix(Graphs, type = "both"))
    } else if (methods::is(Graphs, "list") && methods::is(Graphs[[1]], "igraph")) {
        d <- lapply(Graphs, graphList.to.adjList)
        return(d)
    } else if (methods::is(Graphs, "list") && methods::is(Graphs[[1]], "list")) {
        d <- lapply(Graphs, graphList.to.adjList)
        return(d)
    }
}
