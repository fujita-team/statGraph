#' Graph Spectral Density
#'
#' \code{graph.spectral.density} returns the exact or degree-based spectral density
#' in the interval <\code{from},\code{to}> by using \code{npoints} discretization points.
#'
#' @param Graph the undirected graph (igraph object).
#' If \code{Graph} has the  attribute \code{eigenvalues} containing
#' the eigenvalues of \code{Graph}, such values will be used to
#' compute its spectral density.
#'
#' @param method String that specifies the method to obtain the spectral density. It can
#' take two possible values 'diag' (Default) and 'fast'. If 'diag' is used then
#' the exact spectral density is obtained, otherwise the degree-based spectral
#' density is obtained.
#'
#' @param ... Other relevant parameters to obtain the spectral density such as \code{from}, \code{to},
#' an \code{npoints}.  \code{from}, \code{to} specify  the lower and upper bound of the eigenvalues' support
#' (automatically computed if not given); and \code{npoints} is the number of discretization points (default \code{1024}) of
#' the interval <\code{from},\code{to}>.
#' There are other parameters that depend on the value of the parameter \code{method}:
#' If \code{method='diag'}, then the parameter \code{bandwidth} can be used.
#' This parameter is a string that specifies the criterion to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: 'Silverman' (default), 'Sturges', 'bcv', 'ucv' and 'SJ'.
#' 'bcv' is an abbreviation of biased cross-validation, while 'ucv' means
#' unbiased cross-validation. 'SJ'  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#' Otherwise, if \code{method='fast'}, then the parameter \code{numCores} can be used. This parameter
#' specifies the number of cores (default \code{1}) to use for parallelization.
#'
#'
#' @return  A list with class 'statGraph' containing the following components:
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{info:}}{ a string showing details about the method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' \item{\code{x:}}{ a vector corresponding to the x axis coordinates of the density function.}
#' \item{\code{y:}}{ a vector corresponding to the y axis coordinates of the density function.}
#' \item{\code{from:}}{ a real number corresponding to the smallest value of the x axis.}
#' \item{\code{to:}}{ a real number corresponding to the largest value of the x axis.}
#'
#' @keywords eigenvalue_density
#'
#' @references
#'
#' #' Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
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
#' Newman, M. E. J., Zhang, X., & Nadakuditi, R. R. (2019).
#' Spectra of random networks with arbitrary degrees.
#' Physical Review E, 99(4), 042309.
#'
#' @examples
#' set.seed(1)
#' G <- igraph::sample_smallworld(dim = 1, size = 50, nei = 2, p = 0.2)
#'
#' # Obtain the spectral density
#' density <- graph.spectral.density(Graph = G)
#' density
#'
#' @export
graph.spectral.density <- function(Graph, method = "diag", ...) {
    if (!valid.input(Graph)) {
        stop("The input should be an igraph object!")
    }
    den_fun <- NULL
    if (is.null(Graph$density)) {
        density_parameters <- get.density.parameters(method = method, ...)
        if (method == "diag") {
            den_fun <- graph.diag.spectral.density(Graph = Graph, from = density_parameters$from, to = density_parameters$to, bandwidth = density_parameters$bandwidth, npoints = density_parameters$npoints)
        } else if (method == "fast") {
            den_fun <- graph.fast.spectral.density(Graph = Graph, from = density_parameters$from, to = density_parameters$to, npoints = density_parameters$npoints, numCores = density_parameters$numCores)
        } else {
            stop(paste0(method, " is not a valid method, it only can only be 'diag' or 'fast'."))
        }
    } else {
        den_fun <- Graph$density
    }
    return(den_fun)
}


# Returns the degree-based spectral density in the interval <\code{from},\code{to}> by using npoints discretization points.
graph.fast.spectral.density <- function(Graph, from = NULL, to = NULL, npoints = 1024, numCores = 1) {
    data.name <- deparse(substitute(Graph))
    `%dopar%` <- foreach::`%dopar%`
    # Number of vertices
    n <- igraph::vcount(Graph)
    if (is.null(from)) {
        from <- get.smallest.eigenvalue(Graphs = Graph)
    }
    if (is.null(to)) {
        to <- get.largest.eigenvalue(Graphs = Graph)
    }

    # Discretizise interval <\code{from},\code{to}> in npoints
    bw <- (to - from)/npoints
    x <- seq(from, to, bw)
    y <- rep(0, length(x))
    # Obtain the degree and excess degree distribution
    deg_prob <- c(igraph::degree_distribution(graph = Graph, mode = "all"), 0)  #/vcount(graph)
    k_deg <- seq(1, length(deg_prob)) - 1
    c <- sum(k_deg * deg_prob)
    q_prob <- c()

    for (k in 0:(length(deg_prob) - 1)) {
        aux_q <- (k + 1) * deg_prob[k + 1]/c
        q_prob <- c(q_prob, aux_q)
    }
    # Obtain sorted unique degrees of the graph
    all_k <- c(1:length(q_prob))  #- 1
    valid_idx <- q_prob != 0
    q_prob <- q_prob[valid_idx]
    all_k <- all_k[valid_idx]

    # Obtain the eigenvalue density for each discretized points by using numCores cores.  doMC::registerDoMC(numCores)
    cl <- parallel::makePSOCKcluster(numCores)
    doParallel::registerDoParallel(cl)
    i <- NULL
    y <- foreach::foreach(i = 1:length(x), .combine = c, .export = c("fast.eigenvalue.probability")) %dopar% {
        z <- x[i] + 0.01 * (0 + (0 + (0 + (0+1i))))
        -Im(fast.eigenvalue.probability(deg_prob, q_prob, all_k, z))
    }
    # close cluster
    parallel::stopCluster(cl)
    # return density function as class
    method_info <- "Spectral Density of a Graph"
    info <- "Spectral density obtained with the degree-based method"
    value <- list(method = method_info, info = info, data.name = data.name, x = x, y = y, from = from, to = to)
    class(value) <- "statGraph"
    return(value)
}


# Returns the exact spectral density for a given Graph
graph.diag.spectral.density <- function(Graph, from = NULL, to = NULL, bandwidth = "Silverman", npoints = 1024) {
    data.name <- deparse(substitute(Graph))
    eigenvalues <- graph.eigenvalues(Graph = Graph)
    den_fun <- gaussianDensity(eigenvalues, from, to, bandwidth, npoints)
    ###
    method_info <- "Spectral Density of a Graph"
    info <- "Spectral density obtained with the exact method"
    value <- list(method = method_info, info = info, data.name = data.name, x = den_fun$x, y = den_fun$y, from = den_fun$from, to = den_fun$to)
    class(value) <- "statGraph"
    return(value)
}


# recover parameters relevant for the spectral density relevant for 'method'
get.density.parameters <- function(method = "diag", ...) {
    if (method == "diag")
        return(get.diag.density.parameters(...))
    if (method == "fast")
        return(get.fast.density.parameters(...))
    stop(paste0(method, " is not a valid method, it only can only be 'diag' or 'fast'."))
}

# Returns the probability of an eigenvalue given the degree and excess degree probability.
fast.eigenvalue.probability <- function(deg_prob, q_prob, all_k, z, n_iter = 5000) {
    h_z <- 0 + (0+0i)
    eps <- 1e-07
    all_k_mo <- all_k - 1
    while (n_iter > 0) {
        new_h_z <- sum(q_prob/(1 - all_k_mo * h_z))
        new_h_z <- new_h_z/z^2
        # replaces H_z using the new value found
        if (abs(h_z - new_h_z) < eps) {
            h_z <- new_h_z
            break
        }
        h_z <- new_h_z
        n_iter <- n_iter - 1
    }

    # returns the final result
    count_z <- 0
    for (k in 1:length(deg_prob)) {
        count_z <- count_z + (deg_prob[k])/(1 - (k - 1) * h_z)
    }

    count_z <- (count_z/z)

    return(count_z)
}

# recover parameters relevant for the exact spectral density, or return the default values
get.diag.density.parameters <- function(...) {
    parameters <- list(...)
    valid_parameters <- list()
    param_names <- c("from", "to", "bandwidth", "npoints")
    param_default_vals <- list(NULL, NULL, "Silverman", 1024)
    for (idx in 1:length(param_names)) {
        param_name <- param_names[idx]
        if (is.null(parameters[[param_name]])) {
            valid_parameters[[param_name]] <- param_default_vals[[idx]]
        } else {
            valid_parameters[[param_name]] <- parameters[[param_name]]
        }
    }
    return(valid_parameters)
}

# recover parameters relevant for the degree-based spectral density, or return the default values
get.fast.density.parameters <- function(...) {
    parameters <- list(...)
    valid_parameters <- list()
    param_names <- c("from", "to", "npoints", "numCores")
    param_default_vals <- list(NULL, NULL, 1024, 1)
    for (idx in 1:length(param_names)) {
        param_name <- param_names[idx]
        if (is.null(parameters[[param_name]])) {
            valid_parameters[[param_name]] <- param_default_vals[[idx]]
        } else {
            valid_parameters[[param_name]] <- parameters[[param_name]]
        }
    }
    return(valid_parameters)
}


############################################################################################# SPECTRAL DENSITIES UTILITY FUNCTIONS GIVEN A LIST OF GRAPHS ############

# given a list of graphs that contain their spectral density stored in the 'density' attribute.
get.mean.spectral.density <- function(Graphs) {
    #
    nGraphs <- length(Graphs)
    from <- Graphs[[1]]$density$from
    to <- Graphs[[1]]$density$to
    x <- Graphs[[1]]$density$x
    y <- Reduce(f = "+", Map(f = function(G) {
        G$density$y
    }, Graphs))/nGraphs

    return(list(x = x, y = y, from = from, to = to))
}


# obtains the spectral densities of the graphs in a list.  each spectral density is stored in the 'density' attribute of each graph parameters relevant to the spectral density can be passed
set.list.spectral.density <- function(Graphs, ...) {
    # recover spectral density parameters
    spectral_density_param <- list(...)
    # obtain the extreme eigenvalues
    if (is.null(spectral_density_param$from)) {
        spectral_density_param$from <- get.smallest.eigenvalue(Graphs)
    }

    if (is.null(spectral_density_param$to)) {
        spectral_density_param$to <- get.largest.eigenvalue(Graphs)
    }
    # obtain the spectral density of all graphs
    for (idx in 1:length(Graphs)) {
        spectral_density_param$Graph = Graphs[[idx]]
        Graphs[[idx]] <- do.call(set.spectral.density, spectral_density_param)
    }

    return(Graphs)
}

# computes and sets the spectral density of a graph
set.spectral.density <- function(Graph, ...) {
    Graph$density <- graph.spectral.density(Graph = Graph, ...)
    return(Graph)
}
