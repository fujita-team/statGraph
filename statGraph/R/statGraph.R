#' Graph spectral entropy
#'
#' \code{graph.entropy} returns the spectral entropy of an undirected graph.
#'
#' @param G the undirected graph (igraph type) or its adjacency matrix. The
#' adjacency matrix of an unweighted graph contains only 0s and 1s, while the
#' weighted graph may have nonnegative real values that correspond to the
#' weights of the edges.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @param eigenvalues optional parameter. It contains the eigenvalues of matrix
#' G. Then, if the eigenvalues of matrix G have already been computed, this
#' parameter can be used instead of G. If no value is passed, then the
#' eigenvalues of G will be computed by 'graph.entropy'.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{entropy}{a real number corresponding to the graph spectral entropy.}
#'
#' @keywords spectral_entropy
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
#' G <- igraph::sample_gnp(n=100, p=0.5)
#' entropy <- graph.entropy(G)
#' entropy
#'
#' @export
graph.entropy <- function(G=NULL, bandwidth="Silverman", eigenvalues=NULL) {

  if(methods::is(G,"igraph")){
    A <- as.matrix(igraph::get.adjacency(G))
  }else{
    A <- G
  }
  if (is.null(eigenvalues))
    f <- spectralDensity(A, bandwidth=bandwidth)
  else
    f <- gaussianDensity(eigenvalues, bandwidth=bandwidth)
  if (sum(is.na(f)) > 0)
    return(NA)
  y <- f$y
  i <- which(y != 0)
  y[i] <- y[i]*log(y[i])

  entropy <- -trapezoidSum(f$x, y)

  ###################################################
  method    <- "Spectral entropy of a graph"

  info     <- "Using a Gaussian kernel to estimate the spectral density"
  #info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  value     <- list(method=method, info=info,
                    data.name=data.name, entropy=entropy)
  attr(value, "class") <- "statGraph"

  return(value)
}

#' @export
print.statGraph <- function(x, ...) {
  if(any(names(x) == "method")) cat("\n       ",x$method,"\n\n")

  if(any(names(x) == "info")) cat("-",x$info,"\n\n")

  if(any(names(x) == "data.name")) cat("data:",x$data.name,"\n")

  if(any(names(x) == "entropy")) cat("entropy =", x$entropy, "\n")
  if(any(names(x) == "value")) cat("value =", x$value, "\n")
  if(any(names(x) == "param")) cat("param =", x$param, "\n")
  if(any(names(x) == "KLD")) cat("KLD =", x$KLD, "\n")
  if(any(names(x) == "model")) cat("model:", x$model, "\n")
  if(any(names(x) == "estimates")){
    cat("\nestimates: \n")
    print(x$estimates)
  }
  if(any(names(x) == "values")){
    cat("values:\n")
    print(x$values)
  }

  if(any(names(x) == "cluster")){
    cat("\ncluster:",x$cluster,"\n")
  }
  if(any(names(x) == "parameters")){
    cat("\nparameters:",x$parameters,"\n")
  }

  if(any(names(x) == "x")){
    cat("x:\n")
    print(x$x)
  }
  if(any(names(x) == "y")){
    cat("\n y:\n")
    print(x$y)
  }
  if(any(names(x) == "L1_dist")) cat("L1_dist =", x$L1_dist, "\n")

  cat("\n")
}


#' Graph Information Criterion (GIC)
#'
#' \code{GIC} returns the Kullback-Leibler divergence or L2 distance between an
#' undirected graph and a given graph model.
#'
#' @param G the undirected graph (igraph type) or its adjacency matrix. The
#' adjacency matrix of an unweighted graph contains only 0s and 1s, while the
#' weighted graph may have nonnegative real values that correspond to the
#' weights of the edges.
#'
#' @param model either a list, a string, a function or a matrix describing a
#' graph model:
#'
#' A list that represents the spectral density of a model. It contains the
#' components "x" and "y", which are numeric vectors of the same size. The "x"
#' component contains the points at which the density was computed and the "y"
#' component contains the observed density.
#'
#' A string that indicates one of the following models: "ER" (Erdos-Renyi random
#' graph), "GRG" (geometric random graph), "KR" (k regular random graph), "WS"
#' (Watts-Strogatz model), and "BA" (Barabasi-Albert model). When the argument
#' 'model' is a string, the user must also provides the model parameter by using
#' the argument 'p'.
#'
#' A function that returns a graph (represented by its adjacency matrix)
#' generated by a graph model. It must contain two arguments: the first one
#' corresponds to the graph size and the second to the parameter of the model.
#' The model parameter will be provided by the argument 'p' of the 'GIC'
#' function.
#'
#' A matrix containing the spectrum of the model. Each column contains the
#' eigenvalues of a graph generated by the model. To estimate the spectral
#' density of the model, the method will estimate the density of the values of
#' each column, and then will take the average density.
#'
#' @param p the model parameter. The user must provide a valid parameter if the
#' argument 'model' is a string or a function.
#' For the predefined models ("ER", "GRG", "KR", "WS", and "BA"), the parameter
# 'p' corresponds to:
#
#' the probability to connect a pair of vertices, for the "ER" model
#' (Erdos-Renyi random graph);
#'
#' the radius used to construct the geometric graph in a unit square, for the
#' "GRG" model (geometric random graph);
#'
#' the degree 'k' of a regular graph, for the "KR" model (k regular random
#' graph);
#'
#' the probability to reconnect a vertex, for the "WS" model (Watts-Strogatz
#' model);
#'
#' and the scaling exponent, for the "BA" model (Barabasi-Albert model).
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @param eigenvalues optional parameter. It contains the eigenvalues of matrix
#' G. Then, it can be used when the eigenvalues of G were previously computed.
#' If no value is passed, then the eigenvalues of G will be computed by 'GIC'.
#'
#' @param dist string indicating if you want to use the "KL" (default) or "L2"
#' distances. "KL" means Kullback-Leibler divergence.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{value}{A real number corresponding to the Kullback-Leibler divergence
#' or L2 distance between the observed graph and the graph model..}
#'
#' @keywords graph_information_criterion
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
#' G <- igraph::sample_gnp(n=50, p=0.5)
#'
#' # Using a string to indicate the graph model
#' result1 <- GIC(G, "ER", 0.5)
#' result1
#'
#' # Using a function to describe the graph model
#' # Erdos-Renyi graph
#' model <- function(n, p) {
#'    return(igraph::sample_gnp(n, p))
#' }
#' result2 <- GIC(G, model, 0.5)
#' result2
#' @export
GIC <- function(G, model, p=NULL, bandwidth="Silverman", eigenvalues=NULL,
                dist = "KL") {

  if(methods::is(G,"igraph")){
    A <- as.matrix(igraph::get.adjacency(G))
  }else{
    A <- G
  }
  if (is.null(eigenvalues)){
    eigenvalues <- as.numeric(eigen(A, only.values = TRUE,symmetric=TRUE)$values)
    eigenvalues <- eigenvalues/sqrt(nrow(A))
  }

  if (methods::is(model,"list")) {
    f2 <- model
    f1 <- gaussianDensity(eigenvalues, from=min(f2$x), to=max(f2$x),
                          bandwidth=bandwidth, npoints=1024)
  }else if (methods::is(model,"matrix")) {
    f2 <- nDensities(model, from=min(eigenvalues),
                     to=max(eigenvalues), bandwidth=bandwidth,
                     npoints=1024)
    if (sum(is.na(f2)) > 0)
      return(Inf)
    f2 <- list("x"=f2$x, "y"=rowMeans(f2$densities))
    f1 <- gaussianDensity(eigenvalues, from=min(f2$x), to=max(f2$x),
                          bandwidth=bandwidth, npoints=1024)
  }else {

    fun <- model
    if (methods::is(model,"character")) {
      if (model == "WS"){
        fun <- WSfun(as.integer(sum(A)/(2*ncol(A))))
      }else if(model == "BA"){
        fun <- BAfun(ceiling(mean(rowSums(A))/2))
      }else{
	      fun <- matchFunction(model)
	    }
    }
    f2 <- modelSpectralDensity(fun, ncol(A), p, from=min(eigenvalues),
                               to=max(eigenvalues), bandwidth=bandwidth,
                               npoints=1024)
    if (sum(is.na(f2)) > 0)
      return(Inf)
    f1 <- gaussianDensity(eigenvalues, from=min(f2$x), to=max(f2$x),
                          bandwidth=bandwidth, npoints=1024)
  }
  if (sum(is.na(f1)) > 0)
    return(Inf)
  if (sum(is.na(f2)) > 0)
    return(Inf)
  if(dist == "KL") out <- KL(f1,f2)
  else if(dist == "L2") out <- distance(f1,f2)


  #################################
  method <- "Graph Information criterion"

  info <- "Returns Kullback-Leibler divergence or
  L2 distance between an undirected graph and a given graph model"
  #info <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  output    <- list(method=method, info=info,
                    data.name=data.name, value=out)
  attr(output, "class") <- "statGraph"

  return(output)


  return (out)
}

#' Graph parameter estimator
#'
#' \code{graph.param.estimator} estimates the parameter that best approximates
#' the model to the observed graph according to the Graph Information Criterion
#' (GIC).
#'
#' @param G the undirected graph (igraph type) or its adjacency matrix. The
#' adjacency matrix of an unweighted graph contains only 0s and 1s, while the
#' weighted graph may have nonnegative real values that correspond to the
#' weights of the edges.
#'
#' @param model either a string or a function:
#'
#' A string that indicates one of the following models: "ER" (Erdos-Renyi random
#' graph), "GRG" (geometric random graph), "KR" (k regular random graph), "WS"
#' (Watts-Strogatz model), and "BA" (Barabasi-Albert model).
#'
#' A function that returns a graph (represented by its adjacency matrix)
#' generated by a graph model. It must contain two arguments: the first one
#' corresponds to the graph size and the second to the parameter of the model.
#'
#' @param spectra optional parameter containing the precomputed spectrum of the
#' model. It is a three-dimensional array in which the first dimension
#' corresponds to all parameters that will be explored in the grid, the second
#' dimension has the same size of the given graph, and the third one corresponds
#' to graphs randomly generated by the model. Thus, the position (i,j,k)
#' contains the j-th eigenvalue of the k-th graph generated with the i-th
#' parameter. The attribute 'rownames' of the array corresponds to the
#' parameters converted to string. If spectra is NULL (default), then 'model' is
#' used to generate random graphs and their spectra are computed automatically.
#'
#' @param parameters numeric vector containing the values that that will be
#' considered for the parameter estimation. The 'graph.param.estimator' will
#' return the element of 'parameter' that minimizes the Kullback-Leiber divergence
#' or L1 norm.
#' If the user does not provide the argument 'parameters', and 'model' is an
#' array, then the values considered for the parameter estimation are the
#' rownames converted to a numeric vector. If 'model' is a string, then
#' default values are used for the predefined models ("ER", "GRG", "KR", "WS",
#' and "BA"). The default 'parameter' argument corresponds to a sequence from
#'
#' 0 to 1 with step 'eps' for the "ER" model (Erdos-Renyi random graph), in
#' which the parameter corresponds to the probability to connect a pair of
#' vertices;
#'
#' 0 to sqrt(2) with step 'eps' for the "GRG" model (geometric random graph), in
#' which the parameter corresponds to the radius used to construct the geometric
#' graph in a unit square;
#'
#' 0 to 'n' with step 'n*eps' for the "KR" model (k regular random graph), in
#' which the parameter of the model corresponds to the degree 'k' of a regular
#' graph;
#'
#' 0 to 1 with step 'eps' for the "WS" model (Watts-Strogatz model), in which
#' the parameter corresponds to the probability to reconnect a vertex;
#'
#' and 0 to 3 with step 'eps' for the "BA" model (Barabasi-Albert model), in
#' which the parameter corresponds to the scaling exponent.
#'
#' The default 'parameter' is also used when \code{classic} = FALSE.
#'
#' @param eps precision of the grid (default is 0.01) when 'classic' is TRUE.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @param eigenvalues optional parameter. It contains the eigenvalues of matrix
#' G. Then, it can be used when the eigenvalues of G were previously computed.
#' If no value is passed, then the eigenvalues of G will be computed by
#' 'graph.param.estimator'.
#'
#' @param classic logical. If FALSE parameter is estimated using the fast graph
#' parameter estimator, where this option works better for large graphs with
#' 1000 or more nodes. If TRUE (default) parameter is estimated using grid
#' search.
#'
#' @param npoints Number of points to discretize the interval of the largest and smallest eigenvalue of the
#' graph's adjacency matrix eigenvalues. Only used when \code{classic} = FALSE.
#'
#' @param numCores Number of cores to use for parallelization. Only used when \code{classic} = FALSE.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{param}{the parameter estimate. For the "ER", "GRG", "KR", "WS", and "BA"
#' models, the parameter corresponds to the probability to connect a pair of
#' vertices, the radius used to construct the geometric graph in a unit square,
#' the degree k of a regular graph, the probability to reconnect a vertex, and
#' the scaling exponent, respectively.}
#' \item{KLD}{the Kullback-Leibler divergence between the observed graph and
#' the graph model with the estimated parameter. Only returned if
#' \code{classic} = TRUE.}
#'
#' @keywords parameter_estimation
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
#' G <- igraph::sample_gnp(n=50, p=0.5)
#'
#' # Using a string to indicate the graph model
#' result1 <- graph.param.estimator(G, "ER", eps=0.25)
#' result1
#'
#' \dontrun{
#' # Using a function to describe the graph model
#' # Erdos-Renyi graph
#' set.seed(1)
#' model <- function(n, p) {
#'   return(igraph::sample_gnp(n, p))
#' }
#' result2 <- graph.param.estimator(G, model,  seq(0.2, 0.8, 0.1))
#' result2
#'
#'
#' ### Example giving only the name of the model to use
#' G <- igraph::sample_smallworld(dim = 1, size = 15, nei = 2, p = 0.2)
#'
#' # Obtain the parameter of the WS model
#' result3 <- graph.param.estimator(G, "WS", eps = 1e-1, npoints = 10,
#'                                  numCores = 1)
#' result3
#' }
#'
#' @export
graph.param.estimator <- function(G, model, parameters=NULL, eps=0.01,
                                  bandwidth="Silverman", eigenvalues=NULL,
                                  spectra = NULL, classic = TRUE, npoints = 2000, numCores = 1) {

  if(methods::is(G,"igraph")){
    A <- as.matrix(igraph::get.adjacency(G))
  }else{
    A <- G
  }
  n <- ncol(A)

  if(methods::is(model,"function") && is.null(parameters)) {
    stop("It is necessary to enter the parameters that will be evaluated.")
  }

  if(methods::is(model,"function") && classic == FALSE){
    warning("Changing for classic = TRUE.")
    classic = TRUE
  }

  if (is.null(eigenvalues))
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                      sqrt(nrow(A)))

  #if model is ER, we compute the exact parameter
  if(methods::is(model,"character") && model == "ER"){
    p = sum(A)/(n*(n-1))
    kl <- GIC(A, model, p, bandwidth, eigenvalues=eigenvalues)$value
    out <- list("param" = p, "KLD" = kl)
  }
  #if model is KR, we compute the exact parameter
  else if(methods::is(model,"character") && model == "KR"){
    p = sum(A)/(n)
    kl <- GIC(A, model, p, bandwidth, eigenvalues=eigenvalues)$value
    out <- list("param" = p, "KLD" = kl)
  }
  else if(classic)
  {
    if (!is.null(spectra)) {
      if (is.null(parameters)) parameters <- as.numeric(rownames(spectra))
    }
    else if (is.null(parameters)) {
      if (model == "GRG") parameters <- seq(sqrt(2),0.1, -eps)
      else if (model == "BA") parameters <- seq(3,0.1, -eps)
      else if (model == "KR") parameters <- as.integer(seq(0, 1, eps)*n)
      else parameters <- seq(1, 0.1, -eps)
    }
    pmin <- -1
    klmin <- Inf
    y = c()
    for (p in parameters){
      if (!is.null(spectra)) kl <- GIC(A, spectra[as.character(p),,], p, bandwidth, eigenvalues=eigenvalues)$value
      else kl <- GIC(A, model, p, bandwidth, eigenvalues=eigenvalues)$value
      y = append(y,kl)
      if (kl < klmin) {
        klmin <- kl
        pmin <- p
      }
    }
    out <- list("param"=pmin, "KLD"=klmin)
  }
  else
  {
    out <- fast.graph.param.estimator(G, model, eps = eps,npoints = npoints, numCores = numCores)
  }

  #############################
  method    <- "Graph parameter estimator"

  info     <- "Estimating the parameter that best approximates the model to the observed graph"
  #info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  output     <- list(method=method, info=info,
                     data.name=data.name, param=out$param, KLD=out$KLD)
  attr(output, "class") <- "statGraph"

  return(output)
}

#' Graph model selection
#'
#' \code{graph.model.selection} selects the graph model that best approximates the
#' observed graph according to the Graph Information Criterion (GIC).
#'
#' @param G the undirected graph (igraph type) or its adjacency matrix. The
#' adjacency matrix of an unweighted graph contains only 0s and 1s, while the
#' weighted graph may have nonnegative real values that correspond to the
#' weights of the edges.
#'
#' @param models either a vector of strings, a list of functions or a list of
#' arrays describing graph models:
#'
#' A vector of strings containing some of the following models: "ER" (Erdos-Renyi
#' random graph), "GRG" (geometric random graph), "KR" (k regular random graph),
#' "WS" (Watts-Strogatz model), and "BA" (Barabasi-Albert model).
#'
#' A list of functions. Each function returns a graph (represented by its
#' adjacency matrix) generated by a graph model and has two arguments: the graph
#' size and the model parameter, in this order.
#'
#' A list of arrays. Each elememt of the list is a three-dimensional array
#' containing the precomputed spectrum of each model. Let M be a graph model.
#' For each parameter p considered for M, the array of model M contains the
#' eigenvalues of graphs randomly generated by M with parameter p. The position
#' (i,j,k) of the array contains the j-th eigenvalue of the k-th graph
#' that generated by M with the i-th parameter. The attribute 'rownames' of
#' the array corresponds to the parameters converted to string.
#'
#' If the argument "models" is NULL, then the "ER", "WS", and "BA" models will
#' be considered for the model selection.
#'
#' @param parameters list of numeric vectors. Each vector contains the values
#' that will be considerated for the parameter estimation of the corresponding
#' model.
#' If the user does not provide the argument 'parameters', then default values
#' are used for the predefined models ("ER", "GRG", "KR", "WS", and "BA").
#' The default vector corresponds to a sequence from
#'
#' 0 to 1 with step 'eps' for the "ER" model (Erdos-Renyi random graph), in
#' which the parameter corresponds to the probability to connect a pair of
#' vertices;
#'
#' 0 to sqrt(2) with step 'eps' for the "GRG" model (geometric random graph), in
#' which the parameter corresponds to the radius used to contruct the geometric
#' graph in a unit square;
#'
#' 0 to 'n' with step 'n*eps' for the "KR" model (k regular random graph), in
#' which the parameter of the model corresponds to the degree 'k' of a regular
#' graph;
#'
#' 0 to 1 with step 'eps' for the "WS" model (Watts-Strogatz model), in which
#' the parameter corresponds to the probability to reconnect a vertex;
#'
#' and 0 to 3 with step 'eps' for the "BA" model (Barabasi-Albert model), in
#' which the parameter corresponds to the scaling exponent.
#'
#' @param eps precision of the grid (default is 0.01).
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @param eigenvalues optional parameter. It contains the eigenvalues of matrix
#' G. Then, it can be used when the eigenvalues of G were previously computed.
#' If no value is passed, then the eigenvalues of G will be computed by
#' 'graph.model.selection'.
#'
#' @param ... additional arguments to be used for \code{graph.param.estimator()}
#' function.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{model}{the indice(s) or name(s) of the selected model(s), i. e. the
#' model(s) that minimize(s) the Graph Information Criterion (GIC).}
#' \item{estimates}{a matrix in which each row corresponds to a model, the
#' column "param" corresponds to the parameter estimate, and the column "GIC"
#' corresponds to the Graph Information Criterion (GIC), i. e. the
#' Kullback-Leibler divergence between the observed graph and the model.}
#'
#' @keywords model_selection
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
#'
#' ## Example using an igraph object as input data
#' set.seed(1)
#' G <- igraph::sample_gnp(n=30, p=0.5)
#'
#' # Using strings to indicate the graph models
#' result1 <- graph.model.selection(G, models=c("ER", "WS"), eps=0.5)
#' result1
#'
#'
#' ## Using functions to describe the graph models
#' # Erdos-Renyi graph
#' model1 <- function(n, p) {
#'   return(igraph::sample_gnp(n, p))
#' }
#' # Watts-Strogatz small-world graph
#' model2 <- function(n, pr, K=8) {
#'   return(igraph::sample_smallworld(1, n, K, pr))
#' }
#' parameters <- list(seq(0.01, 0.99, 0.49), seq(0.01, 0.99, 0.49))
#' result2 <- graph.model.selection(G, list(model1, model2), parameters)
#' result2
#'
#' @export
graph.model.selection <- function(G, models=NULL, parameters=NULL, eps=0.01,
                                  bandwidth="Silverman", eigenvalues=NULL, ...) {

  if(methods::is(G,"igraph")){
    A <- as.matrix(igraph::get.adjacency(G))
  }else{
    A <- G
  }
  n <- ncol(A)


  if(methods::is(models,"list") && is.null(parameters)) {
    stop("It is necessary to enter the parameters that will be evaluated.")
  }
  if(is.null(models)){
    models <- c("ER", "WS", "BA")
  }
  results <- matrix(NA, length(models), 2)
  colnames(results) <- c("param", "GIC")


  if (methods::is(models,"character")) {
    rownames(results) <- models
  }


  if(methods::is(models,"list")) {
    if (!is.null(names(parameters)))
      rownames(results) <- names(parameters)
  }
  p <- NULL
  if (is.null(eigenvalues))
    eigenvalues <- (as.numeric(eigen(A, only.values=TRUE,
                                     symmetric=TRUE)$values)/sqrt(nrow(A)))
  for (i in 1:length(models)) {
    if (!is.null(parameters))
      p <- parameters[[i]]

    r <- graph.param.estimator(G, models[[i]], p, eps, bandwidth,
                               eigenvalues=eigenvalues, ...)
    results[i, "param"] <- r$p
    results[i, "GIC"]   <- r$KLD
  }
  m <- which(results[, "GIC"] == min(results[, "GIC"]))
  if (!is.null(rownames(results)))
    m <- rownames(results)[m]

  #############################
  method    <- "Graph model selection"

  info     <- "Selects the graph model that best approximates the observed graph"
  #info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  output     <- list(method=method, info=info,
                     data.name=data.name, model=m, estimates=results)
  attr(output, "class") <- "statGraph"

  return(output)
}

#' Test for the Jensen-Shannon divergence between graphs
#'
#' \code{takahashi.test} tests whether two sets of graphs were generated by the same
#' random graph model.
#' This bootstrap test is based on the Jensen-Shannon (JS) divergence between
#' graphs.
#'
#' Given two lists of graphs, 'G1' and 'G2', 'takahashi.test' tests H0: "JS
#' divergence between 'G1' and 'G2' is 0" against H1: "JS divergence between
#' 'G1' and 'G2' is larger than 0".
#'
#' @param G1 a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param G2 a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param maxBoot integer indicating the number of bootstrap resamplings.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the value of the Jensen-Shannon divergence (JSD) between 'G1' and 'G2'.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @keywords graph_comparison
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
#' G1 <- G2 <- list()
#' for (i in 1:20) {
#'   G1[[i]] <- igraph::sample_gnp(n=50, p=0.5)
#' }
#' for (i in 1:20) {
#'   G2[[i]] <- igraph::sample_gnp(n=50, p=0.51)
#' }
#' result <- takahashi.test(G1, G2, maxBoot=100)
#' result
#'
#' @export
takahashi.test <- function(G1, G2, maxBoot=1000, bandwidth="Silverman") {

  data.name <- paste(deparse(substitute(G1)),"and",deparse(substitute(G2)))

  x <- G1
  y <- G2

  if(methods::is(x,"list") && methods::is(x[[1]],"igraph")){
    x <- f.transform(x)
  }
  if(methods::is(y,"list") && methods::is(y[[1]],"igraph")){
    y <- f.transform(y)
  }

  adjacencyMatrices <- append(x, y)
  labels <- c(rep(0, length(x)), rep(1, length(y)))
  f <- nSpectralDensities(adjacencyMatrices, bandwidth=bandwidth)
  densities <- f$densities
  x <- f$x
  y1 <- rowMeans(densities[, labels==0])
  y2 <- rowMeans(densities[, labels==1])
  n1 <- length(which(labels==0))
  n2 <- length(which(labels==1))
  results <- vector(length=maxBoot)
  ngraphs <- length(adjacencyMatrices)
  result <- JS(list("x"=x, "y"=y1), list("x"=x, "y"=y2))
  for (i in 1:maxBoot) {
    b1 <- sample(1:ngraphs, n1, replace=TRUE)
    b2 <- sample(1:ngraphs, n2, replace=TRUE)
    y1 <- rowMeans(densities[, b1])
    y2 <- rowMeans(densities[, b2])
    results[i] <- JS(list("x"=x, "y"=y1), list("x"=x, "y"=y2))
  }
  pvalue <- (sum(results >= result))/maxBoot


  method           <- "Jensen-Shannon divergence between graphs"
  statistic        <- result
  names(statistic) <- "JSD"
  rval             <- list(statistic=statistic, p.value=pvalue, method=method,
                      data.name=data.name)
  class(rval)      <- "htest"
  return(rval)
}

#' ANOGVA Analysis Of Graph Variability
#'
#' \code{anogva} statistically tests whether two or more sets of graphs are generated
#' by the same random graph model. It is a generalization of the 'takahashi.test'
#' function.
#'
#' @param G a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param labels an array of integers indicating the labels of each graph.
#'
#' @param maxBoot integer indicating the number of bootstrap resamplings.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the statistic of the test.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @keywords analysis_of_graph_variability
#'
#' @examples
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
#'
#' @export
anogva <- function(G, labels, maxBoot=1000, bandwidth="Silverman") {

  data.name <- deparse(substitute(G))

  if(methods::is(G,"list") && methods::is(G[[1]],"igraph")){
    G <- f.transform(G)
  }

  f <- nSpectralDensities(G, bandwidth=bandwidth)
  densidade <- f$densities
  x_axis <- f$x
  band <- length(x_axis)

  media <- matrix(0, max(labels), band)
  mediaAll <- list()
  mediaAll$y <- array(0,band)
  for (j in 1:band) {
    for (i in 1:max(labels)) {
      media[i,j] <- mean(densidade[j, which(labels==i)])
    }
    mediaAll$y[j] <- mean(media[,j])
  }
  mediaAll$x <- x_axis

  distOrig <- 0
  meanGroup <- list()
  meanGroup$x <- x_axis

  for(i in 1:max(labels)) {
    meanGroup$y <- media[i,]
    distOrig <- distOrig + KL(meanGroup, mediaAll)
  }
  distOrig <- distOrig / max(labels)

  ## Permutation test
  distBoot <- array(0, maxBoot)

  for (boot in 1:maxBoot) {
    labelsB <- sample(labels, length(labels), replace=FALSE)

    mediaB <- matrix(0, max(labels), band)

    for (j in 1:band) {
      for (i in 1:max(labels)) {
        mediaB[i,j] <- mean(densidade[j, which(labelsB==i)])
      }
    }

    meanBoot <- list()
    meanBoot$x <- x_axis
    for(i in 1:max(labelsB)) {
      meanBoot$y <- mediaB[i,]
      distBoot[boot] <- distBoot[boot] + KL(meanBoot, mediaAll)
    }
    distBoot[boot] <- distBoot[boot] / max(labels)
  }

  pvalue <- length(which(distBoot >= distOrig)) / (maxBoot+1)

  #htest method
  statistic        <- distOrig
  names(statistic) <- "statistic"
  method           <- "Analysis of Graph Variability"
  rval             <- list(statistic=statistic, p.value=pvalue, method=method,
                           data.name=data.name)
  class(rval)      <- "htest"
  return(rval)
}

#' Semi-Parametric Analysis Of Graph Variability (ANOGVA)
#'
#' \code{sp.anogva} statistically tests whether two or more graphs are
#' generated by the same model and set of parameters.
#'
#' @param G a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param model A string that indicates one of the following models: "ER"
#' (Erdos-Renyi random graph model), "GRG" (geometric random graph model), "WS"
#' (Watts-Strogatz random graph model), and "BA" (Barabasi-Albert random graph
#' model).
#'
#' @param maxBoot integer indicating the number of bootstrap resamples (default
#' is 500).
#'
#' @param spectra optional parameter containing the precomputed spectrum of the
#' model. It is a three-dimensional array in which the first dimension
#' corresponds to all parameters that will be explored in the parameter
#' estimation, the second dimension has the same size of the given graph, and
#' the third one corresponds to graphs randomly generated by the model. Thus,
#' the position (i,j,k) contains the j-th eigenvalue of the k-th graph generated
#' with the i-th parameter. The attribute 'rownames' of the array corresponds to
#' the parameters converted to string. If spectra is NULL (default), then
#' model' is used to generate random graphs and their spectra are computed
#' automatically.
#'
#' @param eps (default is 0.01) precision of the grid when 'classic' = TRUE.
#'
#' @param classic logical. If FALSE parameter is estimated using the fast graph
#' parameter estimator, where this option works better for large graphs with
#' 1000 or more nodes. If TRUE (default) parameter is estimated using grid
#' search.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the F statistic of the test.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{estimates}{a vector containing the estimated parameters for each graph.}
#'
#' @keywords semi_parametric_analysis_of_graph_variability
#'
#' @references
#'
#' Andre Fujita, Eduardo Silva Lira, Suzana de Siqueira Santos, Silvia Yumi
#' Bando, Gabriela Eleuterio Soares, Daniel Yasumasa Takahashi. A
#' semi-parametric statistical test to compare complex networks, Journal of
#' Complex Networks, cnz028, https://doi.org/10.1093/comnet/cnz028
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' _Journal of the Royal Statistical Society series B_, 53, 683-690.
#' http://www.jstor.org/stable/2345597.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' model <- "ER"
#' G <- list()
#'
#' # Under H0
#' G[[1]] <- igraph::sample_gnp(50, 0.5)
#' G[[2]] <- igraph::sample_gnp(50, 0.5)
#' G[[3]] <- igraph::sample_gnp(50, 0.5)
#' result1 <- sp.anogva(G, model, maxBoot = 300)
#' result1
#'
#' # Under H1
#' G[[1]] <- igraph::sample_gnp(50, 0.5)
#' G[[2]] <- igraph::sample_gnp(50, 0.55)
#' G[[3]] <- igraph::sample_gnp(50, 0.5)
#' result2 <- sp.anogva(G, model, maxBoot = 300)
#' result2
#' }
#'
#' @export
sp.anogva <- function(G, model, maxBoot=500, spectra = NULL, eps = 0.01,
                      classic = TRUE, bandwidth = "Silverman") {

  data.name <- deparse(substitute(G))
  graph     <- G

  if(methods::is(graph,"list") && methods::is(graph[[1]],"igraph")){
    graph <- f.transform(graph)
  }

  g <- length(graph)
  p.hat <- list()
  for(l in 1:g) {
    n <- as.character(ncol(graph[[l]]))

    if(methods::is(spectra,"list"))
      p.hat[[l]] <- graph.param.estimator(graph[[l]], model, eps=eps,
                                          bandwidth = bandwidth,
                                          spectra = spectra[[n]],
                                          classic = classic)$p
    else p.hat[[l]] <- graph.param.estimator(graph[[l]], model, eps=eps,
                                             bandwidth = bandwidth,
                                             spectra = spectra,
                                             classic = classic)$p
  }
  p.boot <- matrix(0, length(graph), maxBoot)
  for (boot in 1:maxBoot) {
    g.boot <- list()
    for (l in 1:g) {
      n <- as.character(ncol(graph[[l]]))
      if (model == "ER") {
        g.boot[[l]] <- ER(ncol(graph[[l]]), p.hat[[l]])
      }
      else if (model == "GRG") {
        g.boot[[l]] <- GRG(ncol(graph[[l]]), p.hat[[l]])
      }
      else if (model == "WS") {
        g.boot[[l]] <- WS(ncol(graph[[l]]),p.hat[[l]],2)
      }
      else if (model == "BA") {
        g.boot[[l]] <- BA(ncol(graph[[l]]), p.hat[[l]])
      }
      if(methods::is(spectra,"list"))
        p.boot[l,boot] <- graph.param.estimator(g.boot[[l]], model, eps=eps,
                                                bandwidth = bandwidth,
                                                spectra = spectra[[n]],
                                                classic = classic)$p
      else p.boot[l,boot] <- graph.param.estimator(g.boot[[l]], model, eps=eps,
                                                   bandwidth = bandwidth,
                                                   spectra = spectra,
                                                   classic = classic)$p
    }
  }
  var.boot <- array(0, g)
  for(l in 1:g) {
    var.boot[l] <- var(p.boot[l,])
  }
  SSres <- 0
  SStr <- 0
  m <- mean(as.numeric(as.array(p.hat)))
  for(l in 1:g) {
    SSres <- SSres + (maxBoot-1) * var.boot[l]
    SStr <- SStr + (p.hat[[l]] - m)^2
  }
  F_ <- (SStr / (g-1)) / ( SSres / (g*maxBoot - g))
  p <- pf(F_, df1=(g-1), df2=(g*maxBoot - g), lower.tail=FALSE)
  res <- list()
  res$parameters <- unlist(p.hat) #as.array(p.hat)
  res$F.value <- F_
  res$p.value <- p


  #htest method
  statistic         <- res$F.value
  names(statistic)  <- "F.value"
  #parameter        <- res$parameters
  #names(parameter) <- paste("parameter",1:length(parameter),sep='')
  estimate          <- res$parameters
  names(estimate)   <- paste("parameter",1:length(estimate),sep='')
  method            <- "Semi-Parametric Analysis Of Graph Variability"
  rval              <- list(statistic=statistic,
                            p.value=res$p.value, method=method,
                            data.name=data.name, estimate=estimate)
  class(rval)       <- "htest"
  return(rval)
}

#' Test for Association / Correlation Between Paired Samples of Graphs
#'
#' \code{graph.cor.test} tests for association between paired samples of graphs,
#' using Spearman's rho correlation coefficient.
#'
#' @param G1 a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param G2 a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{estimates}{the estimated measure of association 'rho'.}
#'
#' @keywords correlation_coefficient
#'
#' @references
#' Fujita, A., Takahashi, D. Y., Balardin, J. B., Vidal, M. C. and Sato, J. R.
#' (2017) Correlation between graphs with an application to brain network
#' analysis. _Computational Statistics & Data Analysis_ *109*, 76-92.
#'
#' @examples
#' set.seed(1)
#' G1 <- G2 <- list()
#'
#' p <- MASS::mvrnorm(50, mu=c(0,0), Sigma=matrix(c(1, 0.5, 0.5, 1), 2, 2))
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
#' @import MASS
#' @export
graph.cor.test <- function(G1, G2) {

  data.name <- paste(deparse(substitute(G1)),"and",deparse(substitute(G2)))

  if(methods::is(G1,"list") && methods::is(G1[[1]],"igraph")){
    G1 <- f.transform(G1)
  }
  if(methods::is(G2,"list") && methods::is(G2[[1]],"igraph")){
    G2 <- f.transform(G2)
  }

  G1.radius <- array(0, length(G1))
  G2.radius <- array(0, length(G2))

  for (i in 1:length(G1)) {
    G1.radius[i] <- eigen(G1[[i]], only.values=TRUE, symmetric=TRUE)$values[1]
    G2.radius[i] <- eigen(G2[[i]], only.values=TRUE, symmetric=TRUE)$values[1]
  }


  res <- cor.test(G1.radius, G2.radius, method="spearman")
  #return(cor.test(G1.radius, G2.radius, method="spearman"))

  #htest method
  statistic         <- res$statistic
  names(statistic)  <- "statistic"
  #parameter        <- res$parameters
  #names(parameter) <- paste("parameter",1:length(parameter),sep='')
  estimate          <- res$estimate
  names(estimate)   <- "rho"
  method            <- "Association between paired samples of graphs, using Spearman's rho correlation coefficient"
  rval              <- list(statistic=statistic,
                            p.value=res$p.value, method=method,
                            data.name=data.name, estimate=estimate)
  class(rval)       <- "htest"
  return(rval)
}

#' Auto Correlation Function Estimation for Graphs
#'
#' The function \code{graph.acf} computes estimates of the autocorrelation
#' function for graphs.
#'
#' @param G a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param plot logical. If TRUE (default) the graph.acf is plotted.
#'
#' @return An object of class acf.
#'
#' @keywords autocorrelation
#'
#' @references
#' Fujita, A., Takahashi, D. Y., Balardin, J. B., Vidal, M. C. and Sato, J. R.
#' (2017) Correlation between graphs with an application to brain network
#' analysis. _Computational Statistics & Data Analysis_ *109*, 76-92.
#'
#' @examples
#' set.seed(1)
#' G <- list()
#' p <- array(0, 100)
#' p[1:3] <- rnorm(3)
#' for (t in 4:100) {
#'   p[t] <- 0.5*p[t-3] + rnorm(1)
#' }
#' ma <- max(p)
#' mi <- min(p)
#' p <- (p - mi)/(ma-mi)
#' for (t in 1:100) {
#'   G[[t]] <- igraph::sample_gnp(100, p[t])
#' }
#' graph.acf(G, plot=TRUE)
#'
#' @import stats
#' @export
graph.acf <- function(G, plot=TRUE) {

  if(methods::is(G,"list") && methods::is(G[[1]],"igraph")){
    G <- f.transform(G)
  }
  G.radius <- array(0, length(G))
  for (t in 1:length(G)) {
    G.radius[t] <- eigen(G[[t]], only.values=TRUE, symmetric=TRUE)$values[1]
  }
  res <- acf(G.radius, plot=plot)
  return(res)
}

#' Hierarchical cluster analysis on a list of graphs.
#'
#' Given a list of graphs, \code{graph.hclust} builds a hierarchy of clusters
#' according to the Jensen-Shannon divergence between graphs.
#'
#' @param G a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param k the number of clusters.
#'
#' @param method the agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of '"ward.D"', '"ward.D2"', '"single"',
#' '"complete"', '"average"' (= UPGMA), '"mcquitty"' (= WPGMA), '"median"'
#' (= WPGMC) or '"centroid"' (= UPGMC).
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @return A list containing:
#' \item{hclust}{an object of class *hclust* which describes the tree produced
#' by the clustering process.}
#' \item{cluster}{the clustering labels for each graph.}
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
#'
#' @export
graph.hclust <- function(G, k, method="complete", bandwidth="Silverman") {
  x <- G

  if(methods::is(x,"list") && methods::is(x[[1]],"igraph")){
    x <- f.transform(x)
  }
  f <- nSpectralDensities(x, bandwidth=bandwidth)

  d <- matrix(0, length(x), length(x))
  for (i in 1:(length(x)-1)) {
    f1 <- list("x"=f$x, "y"=f$densities[,i])
    for (j in (i+1) : length(x)) {
      f2 <- list("x"=f$x, "y"=f$densities[,j])
      d[i,j] <- d[j,i] <- sqrt(JS(f1, f2))
    }
  }
  tmp <- hclust(as.dist(d), method=method)

  res <- list()
  res$hclust <- tmp
  res$cluster <- cutree(tmp, k)

  return(res)
}

#' Multidimensional scaling of graphs
#'
#' \code{graph.mult.scaling} performs multidimensional scaling of graphs. It
#' takes the Jensen-Shannon divergence between graphs (JS) and uses the
#' 'cmdscale' function from the 'stats' package to obtain a set of points such
#' that the distances between the points are similar to JS.
#'
#' @param G a list of undirected graphs (igraph type) or their adjacency
#' matrices. The adjacency matrix of an unweighted graph contains only 0s and
#' 1s, while the weighted graph may have nonnegative real values that correspond
#' to the weights of the edges.
#'
#' @param plot logical. If TRUE (default) the points chosen to represent the
#' Jensen-Shannon divergence between graphs are plotted.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#' @param type what type of plot should be drawn. The defaut value is '"n"',
#' which indicates that the points will not be plotted (i. e. only the labels
#' of the graphs will be plotted).
#'
#' @param main title of the plot (default value is "").
#'
#' @param ... additional plotting parameters. See 'plot' function from the
#' 'graphics' package for the complete list.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{values}{A matrix in which each column corresponds to a coordinate and each
#' row corresponds to a graph (point). Then, each row gives the coordinates of
#' the points chosen to represent the Jensen-Shannon divergence between graphs.}
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
#' @import graphics
#' @export
graph.mult.scaling <- function(G, plot=TRUE, bandwidth="Silverman", type="n",
                               main="", ...) {
  x <- G

  if(methods::is(x,"list") && methods::is(x[[1]],"igraph")){
    x <- f.transform(x)
  }

  f <- nSpectralDensities(x, bandwidth=bandwidth)

  d <- matrix(0, length(x), length(x))
  for (i in 1:(length(x)-1)) {
    f1 <- list("x"=f$x, "y"=f$densities[,i])
    for (j in (i+1) : length(x)) {
      f2 <- list("x"=f$x, "y"=f$densities[,j])
      d[i,j] <- d[j,i] <- sqrt(JS(f1, f2))
    }
  }

  if (is.null(names(x)))
    names <- as.character(1:length(x))
  else
    names <- names(x)
  colnames(d) <- rownames(d) <- names
  fit <- cmdscale(as.dist(d), k=2)

  x <- fit[,1]
  y <- fit[,2]
  names(x) <- names
  names(y) <- names
  if (plot) {
    plot(x, y, xlab="x", ylab="y", main=main, type=type, ...)
    text(x, y, labels=names, cex=1)
  }

  ######################

  method    <- "Multidimensional scaling of graphs"

  #info     <- "Taking the Jensen-Shannon divergence between graphs (JS) and
  #using the 'cmdscale' function from the 'stats' package to obtain a set of
  #points such that the distances between the points are similar to JS"
  info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  output     <- list(method=method, info=info,
                     values=fit)
  attr(output, "class") <- "statGraph"

  return(output)
}

#' Tang hypothesis testing for random graphs.
#'
#' Given two independent finite-dimensional random dot product graphs,
#' \code{tang.test} tests if they have generating latent positions that are drawn
#' from the same distribution.
#'
#' @param G1 the first undirected graph to be compared. Must be an igraph
#' object.
#'
#' @param G2 the second undirected graph to be compared. Must be an igraph
#' object.
#'
#' @param dim dimension of the adjacency spectral embedding.
#'
#' @param sigma a real value indicating the kernel bandwidth. If NULL (default)
#' the bandwidth is calculated by the method.
#'
#' @param maxBoot integer indicating the number of bootstrap resamples
#' (default is 200).
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the T-value of the test.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @references
#' Tang, Minh, et al. "A nonparametric two-sample hypothesis testing problem for
#' random graphs." Bernoulli 23.3 (2017): 1599-1630.
#'
#' Tang, Minh, et al. "A semiparametric two-sample hypothesis testing problem
#' for random graphs." Journal of Computational and Graphical Statistics 26.2
#' (2017): 344-354.
#'
#' @examples
#' set.seed(42)
#'
#' ## test under H0
#' lpvs <- matrix(rnorm(200), 20, 10)
#' lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
#' G1 <- igraph::sample_dot_product(lpvs)
#' G2 <- igraph::sample_dot_product(lpvs)
#' D1 <- tang.test(G1, G2, 5)
#' D1
#'
#' ## test under H1
#' lpvs2 <- matrix(pnorm(200), 20, 10)
#' lpvs2 <- apply(lpvs2, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
#' G2 <- suppressWarnings(igraph::sample_dot_product(lpvs2))
#' D2 <- tang.test(G1, G2, 5)
#' D2
#'
#' @export
tang.test <- function(G1, G2, dim, sigma = NULL, maxBoot=200){

  data.name <- paste(deparse(substitute(G1)),"and",deparse(substitute(G2)))

  t.validateInput(G1, G2, dim, maxBoot)
  Xhat1 = t.embed.graph(G1, dim)
  Xhat2 = t.embed.graph(G2, dim)
  if(is.null(sigma)){
    sigma = t.get.sigma(Xhat1, Xhat2)
  }
  test_stat = t.test.stat(Xhat1, Xhat2, sigma)
  test_distribution = t.sampling.distribution(G1, dim, maxBoot)
  #test_distribution2 = sampling.distribution(G2, dim, maxBoot)
  p_val = t.p_value(test_stat, test_distribution)


  #htest method
  statistic         <- test_stat
  names(statistic)  <- "T"
  method            <- "Tang hypothesis testing for random graphs"
  rval              <- list(statistic=statistic, p.value=p_val, method=method,
                            data.name=data.name)
  class(rval)       <- "htest"
  return(rval)
}

#' Ghoshdastidar hypothesis testing for large random graphs.
#'
#' Given two lists of graphs generated by the inhomogeneous random graph model,
#' \code{ghoshdastidar.test} tests if they were generated by the same parameters.
#'
#' @param G1 the first list of undirected graphs to be compared. Must be a list
#' of matrices or igraph objects.
#'
#' @param G2 the second list of undirected graphs to be compared. Must be a list
#' of matrices or igraph objects.
#'
#' @param maxBoot integer indicating the number of bootstrap resamples (default
#' is 300).
#'
#' @param two.sample logical. If TRUE the sets contain only one graph each. If
#' FALSE the sets contain more than one graph each (default is FALSE).
#'
#
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the T-value of the test.}
#' \item{p.value}{the p-value of the test (only returned when the parameter
#' 'two.sample' is FALSE).}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @references
#' Ghoshdastidar, Debarghya, et al. "Two-sample tests for large random graphs
#' using network statistics". arXiv preprint arXiv:1705.06168 (2017).
#'
#' Ghoshdastidar, Debarghya, et al. "Two-sample hypothesis testing for
#' inhomogeneous random graphs". arXiv preprint, arXiv:1707.00833 (2017).
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#'
#' ## test for sets with more than one graph each under H0
#' G1 <- G2 <- list()
#' for(i in 1:10){
#'   G1[[i]] <- as.matrix(igraph::get.adjacency(igraph::sample_gnp(50,0.6)))
#'   G2[[i]] <- as.matrix(igraph::get.adjacency(igraph::sample_gnp(50,0.6)))
#' }
#' D1 <- ghoshdastidar.test(G1, G2)
#' D1
#'
#' ## test for sets with more than one graph each under H1
#' G1 <- G2 <- list()
#' for(i in 1:10){
#'   G1[[i]] <- as.matrix(igraph::get.adjacency(igraph::sample_gnp(50,0.6)))
#'   G2[[i]] <- as.matrix(igraph::get.adjacency(igraph::sample_gnp(50,0.7)))
#' }
#' D2 <- ghoshdastidar.test(G1, G2)
#' D2
#'
#' ## test for sets with only one graph each under H0
#' G1 <- G2 <- list()
#' G1[[1]] <- igraph::sample_gnp(300, 0.6)
#' G2[[1]] <- igraph::sample_gnp(300, 0.6)
#' D3 <- ghoshdastidar.test(G1, G2, two.sample= TRUE)
#' D3
#'
#' ## test for sets with only one graph each under H1
#' G1 <- G2 <- list()
#' G1[[1]] <- igraph::sample_gnp(300, 0.6)
#' G2[[1]] <- igraph::sample_gnp(300, 0.7)
#' D4 <- ghoshdastidar.test(G1, G2, two.sample= TRUE)
#' D4
#' }
#'
#' @export
ghoshdastidar.test <- function(G1, G2, maxBoot = 300, two.sample = FALSE)
{

  data.name <- paste(deparse(substitute(G1)),"and",deparse(substitute(G2)))

  if(methods::is(G1,"list") && methods::is(G1[[1]],"igraph")){
    G1 <- g.transform(G1)
  }
  if(methods::is(G2,"list") && methods::is(G2[[1]],"igraph")){
    G2 <- g.transform(G2)
    }
  if(methods::is(G1,"igraph") && methods::is(G2,"igraph")){
    G1 <- g.transform(G1)
    G2 <- g.transform(G2)
  }


  if(!methods::is(G1,"list") || !methods::is(G1[[1]],"matrix")) stop("G1 must be a list of matrices or igraph objects.")
  if(!methods::is(G2,"list") || !methods::is(G2[[1]],"matrix")) stop("G2 must be a list of matrices or igraph objects.")

  D <- g.test(G1, G2)


  statistic         <- D
  names(statistic)  <- "T"
  method  <- "Ghoshdastidar hypothesis testing for large random graphs"

  if( !two.sample ){
    test_distribution <- g.sampling.distribution(G1,G2,maxBoot)
    p_val <- mean(test_distribution >= D)
    #out = list(T = D, p.value = p_val)

    rval <- list(statistic=statistic, p.value=p_val,
                 method=method, data.name=data.name)
  }
  else{
    #out = list(T = D)

    rval <- list(statistic=statistic,
                 method=method, data.name=data.name)
  }

  class(rval) <- "htest"
  return(rval)

}

#' Cerqueira et al.’s test
#'
#' Given two identically independently distributed (idd) samples of graphs G1 and
#' G2, the test verifies if they have the same distribution by calculating the
#' mean distance D from G1 to G2. The test rejects the null hypothesis if D is
#' greater than the (1-alpha)-quantile of the distribution of the test under the
#' null hypothesis.
#'
#' @param G1 the first iid sample of graphs to be compared. Must be a list of
#' igraph objects.
#'
#' @param G2 the second iid sample of graphs to be compared. Must be a list of
#' igraph objects.
#'
#' @param maxBoot integer indicating the number of bootstrap resamples (default
#' is 300).
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the W-value of the test.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @references
#' Andressa Cerqueira, Daniel Fraiman, Claudia D. Vargas and Florencia Leonardi.
#' "A test of hypotheses for random graph distributions built from EEG data",
#' https://ieeexplore.ieee.org/document/7862892
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#'
#' ## test under H0
#' G1 <- G2 <- list()
#' for(i in 1:10){
#'   G1[[i]] <- igraph::sample_gnp(50,0.5)
#'   G2[[i]] <- igraph::sample_gnp(50,0.5)
#' }
#' k1 <- cerqueira.test(G1, G2)
#' k1
#'
#' ## test under H1
#' G1 <- G2 <- list()
#' for(i in 1:10){
#'   G1[[i]] <- igraph::sample_gnp(50,0.5)
#'   G2[[i]] <- igraph::sample_gnp(50,0.6)
#' }
#' k2 <- cerqueira.test(G1, G2)
#' k2
#' }
#'
#' @export
cerqueira.test <- function(G1, G2, maxBoot = 300){

  data.name <- paste(deparse(substitute(G1)),"and",deparse(substitute(G2)))

  if(methods::is(G1,"list") && methods::is(G1[[1]],"igraph")){
    G1 <- c.transform(G1)
  }
  else{
    stop("Parameter G1 must be a list of igraph objects.")
  }

  if(methods::is(G2,"list") && methods::is(G2[[1]],"igraph")){
    G2 <- c.transform(G2)
  }
  else{
    stop("Parameter G2 must be a list of igraph objects.")
  }

  D <- c.test(G1,G2)

  test_distribution <- c.sampling.distribution(G1,G2,maxBoot)

  p_val <- mean(test_distribution >= D)

  #htest method
  statistic         <- D
  names(statistic)  <- "W"
  method            <- "Verify if two samples of random graphs were originated from the same probability distribution."
  rval              <- list(statistic=statistic, p.value=p_val, method=method,
                            data.name=data.name)
  class(rval)       <- "htest"
  return(rval)
}


#' Fraiman’s test
#'
#' Given a list of graphs, the test verifies if all the subpopulations have the
#' same mean network, under the alternative that at least one subpopulation has
#' a different mean network.
#'
#' @param G the undirected graphs to be compared. Must be a list of lists of
#' igraph objects or a list of lists of adjacency matrices.
#'
#' @param maxBoot integer indicating the number of bootstrap resamples
#' (default is 300).
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the T-value of the test.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a string indicating the used method.}
#' \item{data.name}{a string with the data's name(s).}
#'
#' @references
#' Fraiman, Daniel, and Ricardo Fraiman. "An ANOVA approach for statistical
#' comparisons of brain networks",
#' https://www.nature.com/articles/s41598-018-23152-5
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#'
#' ## test under H0
#' a <- b <- G <- list()
#' for(i in 1:10){
#'   a[[i]] <- igraph::sample_gnp(50,0.5)
#'   b[[i]] <- igraph::sample_gnp(50,0.5)
#' }
#' G <- list(a,b)
#' k1 <- fraiman.test(G)
#' k1
#'
#' ## test under H1
#' a <- b <- G <- list()
#' for(i in 1:10){
#'   a[[i]] <- igraph::sample_gnp(50,0.5)
#'   b[[i]] <- igraph::sample_gnp(50,0.6)
#' }
#' G <- list(a,b)
#' k2 <- fraiman.test(G)
#' k2
#' }
#'
#' @export
fraiman.test <- function(G, maxBoot = 300){

  data.name <- deparse(substitute(G))

  # transform and verify input
  if(methods::is(G,"list") && methods::is(G[[1]],"list") && methods::is(G[[1]][[1]],"igraph")){
    G <- f.transform(G)
  }
  if(!methods::is(G,"list") || !methods::is(G[[1]],"list") || !methods::is(G[[1]][[1]],"matrix"))
    stop(paste("You must pass a list of lists of igraphs or a list of lists of",
               "matrices."))

  D <- f.test(G)

  test_distribution <- f.sampling.distribution(G, maxBoot)

  ## modification made in april 1st, 2019:
  p_val <- mean(test_distribution <= D)

  ## original was:
  # p_val <- mean(abs(test_distribution) >= abs(D))

  #htest method
  statistic         <- D
  names(statistic)  <- "T"
  method            <- "Test for network differences between groups with an analysis of variance test (ANOVA)"
  rval              <- list(statistic=statistic, p.value=p_val, method=method,
                            data.name=data.name)
  class(rval)       <- "htest"
  return(rval)
}



################################################################################
## Auxiliary functions
################################################################################

# Returns the Jensen-Shannon divergence between two densities
JS <- function(f1, f2) {
  fm <- f1
  fm$y <- (f1$y + f2$y)/2
  return((KL(f1, fm) + KL(f2, fm))/2)
}

# Returns the Kullback-Leibler divergence between two densities
KL <- function(f1, f2) {

  y <- f1$y
  #n <- length(y)
  #for (i in 1:n) {
  #  if (y[i] != 0 && f2$y[i] == 0){
  #    return (Inf)
  #  }
  #  if (y[i] != 0)
  #    y[i] <- y[i]*log(y[i]/f2$y[i])
  #}
  y <- ifelse(y != 0, y*log(y/f2$y), y)

  return (trapezoidSum(f1$x, y))
}

# Given a partition x[1]...x[n] and y[i] = f(x[i]), returns the trapezoid sum
# approximation for int_{x[1]}^{x[n]}{f(x)dx}
trapezoidSum <- function (x, y) {
  n <- length(x)
  delta <- (x[2] - x[1])
  area <- sum(y[2:(n-1)])
  area <- (area + (y[1] + y[n])/2)*delta
  return(area)
}



# Returns the kernel bandwidth for a sample x based on Sturge's criterion
kernelBandwidth <- function(x) {
  n <- length(x)
  nbins <- ceiling(log2(n) + 1)
  return(abs(max(x) - min(x))/nbins)
}


# Returns the density function for a sample x at n points in the interval [from, to]
gaussianDensity <- function(x, from=NULL, to=NULL, bandwidth="Silverman",
                            npoints=1024) {
  if (bandwidth == "Sturges"){
    bw <- kernelBandwidth(x)
  }else if (bandwidth == "Silverman"){
    bw <- bw.nrd0(x)
  }else if (bandwidth == "bcv"){
    bw <- suppressWarnings(bw.bcv(x))
  }else if (bandwidth == "ucv"){
    bw <- suppressWarnings(bw.ucv(x))
  }else if (bandwidth == "SJ"){
    bw <- "SJ"
  }else{
    stop("Please, choose a valid bandwidth.")
  }
  if (bw == 0){
    stop("bw cannot be zero.")
  }
  if (is.null(from) || is.null(to)){
    f <- density(x, bw=bw, n=npoints)
  }else{
    f <- density(x, bw=bw, from=from, to=to, n=npoints)
  }

  f$y = f$y + 1e-12 # we do not want the area to be zero, so we add a very small number
  area <- trapezoidSum(f$x, f$y)
  return(list("x"=f$x, "y"=f$y/area))
}

# Returns the spectral density for a given adjacency matrix A
spectralDensity <- function(A, from=NULL, to=NULL, bandwidth="Silverman",
                            npoints=1024) {
  eigenvalues <- as.numeric(eigen(A, only.values=TRUE, symmetric=TRUE)$values)
  eigenvalues <- eigenvalues/sqrt(nrow(A))
  return(gaussianDensity(eigenvalues, from, to, bandwidth, npoints))
}

nDensities <- function (spectra, from=NULL, to=NULL,
                        bandwidth="Silverman", npoints=1024) {
  ngraphs <- ncol(spectra)
  densities <- matrix(NA, npoints, ngraphs)
  minimum <- min(spectra)
  maximum <- max(spectra)
  if (!is.null(from) && !is.null(to)) {
    minimum <- min(minimum, from)
    maximum <- max(maximum, to)
  }
  for (i in 1:ngraphs) {
    f <- gaussianDensity(spectra[,i], bandwidth=bandwidth,
                         from=minimum, to=maximum,
                         npoints=npoints)
    if (sum(is.na(f)) > 0) {
      return(NA)
    }
    else {
      densities[,i] <- f$y
      x <- f$x
    }
  }
  return(list("x"=x, "densities"=densities))
}


# Returns the spectral densities for given adjacency matrices A1 and A2 at the
# same points
spectralDensities <- function(A1, A2, bandwidth="Silverman",
                              npoints=1024) {
  n1 <- nrow(A1)
  n2 <- nrow(A2)
  e1 <- (as.numeric(eigen(A1, only.values = TRUE,
                          symmetric=TRUE)$values)/sqrt(n1))
  e2 <- (as.numeric(eigen(A2, only.values = TRUE,
                          symmetric=TRUE)$values)/sqrt(n2))
  #b1 <- kernelBandwidth(e1)
  #b2 <- kernelBandwidth(e2)
  #from <- min(min(e1) - 3*b1, min(e2) - 3*b2)
  #to <- max(max(e1) + 3*b1, max(e2) + 3*b2)
  from <- min(e1, e2)
  to <- max(e1, e2)
  f1 <- gaussianDensity(e1, from=from, to=to, bandwidth=bandwidth,
                        npoints=npoints)
  f2 <- gaussianDensity(e2, from=from, to=to, bandwidth=bandwidth,
                        npoints=npoints)
  if (sum(is.na(f1)) > 0 || sum(is.na(f2)) > 0)
    return(NA)
  return(list("f1"=f1, "f2"=f2))
}

# Returns the spectral densities for a list of adjacency matrices at the
# same points
nSpectralDensities <- function (adjacencyMatrices, from=NULL, to=NULL,
                                bandwidth="Silverman") {
  npoints <- 1024
  ngraphs <- length(adjacencyMatrices)
  ns <- unlist(lapply(adjacencyMatrices, ncol))
  #n <- ncol(adjacencyMatrices[[1]])
  spectra <- matrix(NA, max(ns), ngraphs)
  for (i in 1:ngraphs) {
    A <- adjacencyMatrices[[i]]
    n <- ncol(A)
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE,
                                     symmetric=TRUE)$values)/sqrt(n))
    spectra[1:n,i] <- eigenvalues
  }
  densities <- matrix(NA, npoints, ngraphs)
  minimum <- min(spectra, na.rm=TRUE)
  maximum <- max(spectra, na.rm=TRUE)
  if (!is.null(from) && !is.null(to)) {
    minimum <- min(minimum, from)
    maximum <- max(maximum, to)
  }
  for (i in 1:ngraphs) {
    n <- ns[i]
    f <- gaussianDensity(spectra[1:n,i], bandwidth=bandwidth,
                         from=minimum, to=maximum,
                         npoints=npoints)

    densities[,i] <- f$y
    x <- f$x
  }
  return(list("x"=x, "densities"=densities))
}

# Estimates the spectral density of a graph model
modelSpectralDensity <- function(fun, n, p, ngraphs=100, from=NULL, to=NULL,
                                 bandwidth="Silverman", npoints=1024) {
  spectra <- matrix(NA, n, ngraphs)
  for (i in 1:ngraphs) {
    A <- fun(n, p)


    if(methods::is(A,"igraph")){
      A <- as.matrix(igraph::get.adjacency(A))
    }
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE,
                                     symmetric=TRUE)$values)/sqrt(nrow(A)))
    spectra[,i] <- eigenvalues
  }
  densities <- matrix(NA, npoints, ngraphs)
  minimum <- min(spectra)
  maximum <- max(spectra)
  if (!is.null(from) && !is.null(to)) {
    minimum <- min(minimum, from)
    maximum <- max(maximum, to)
  }
  for (i in 1:ngraphs) {
    f <- gaussianDensity(spectra[,i], from=minimum, to=maximum,
                         bandwidth=bandwidth, npoints=npoints)

    densities[,i] <- f$y
    x <- f$x
  }
  return(list("x" = x, "y" = rowMeans(densities)))
}

modelSpectra <- function(model, n, p, ngraphs=50) {
  fun <- model
  if (is.character(model)) {
    if (model == "WS")
      fun <- WSfun(2)
    else fun <- matchFunction(model)
  }
  spectra <- matrix(NA, n, ngraphs)
  for (i in 1:ngraphs) {
    A <- fun(n, p)
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE,
                                     symmetric=TRUE)$values)/sqrt(nrow(A)))
    spectra[,i] <- eigenvalues
  }
  return(spectra)
}

# Extract a function specified by name
matchFunction <- function(name) {
  return(match.fun(name))
}


## Auxiliary for Cerqueira method. Test distribution under the null hypothesis
c.sampling.distribution <- function(g, gp, maxBoot = 300)
{

  m <- nrow(g)+nrow(gp)
  test_distribution = c()
  for (i_per in 1:maxBoot){
    total <- rbind(g, gp)
    ind <- sample(1:m, floor(m/2), replace=F)
    xa <- total[ind,]
    ya <- total[-ind,]
    test_distribution[i_per] <- c.test(xa,ya)
  }
  return(sort(test_distribution))
}

## Auxiliary for Cerqueira method. Fix input format.
c.transform <- function(g, n = igraph::gorder(g[[1]]))
{
  x <- matrix(0, length(g), n*(n-1)/2)
  i <- 1
  for(gr in g){
    aux <- as.matrix(igraph::get.adjacency(gr))
    x[i,] <- aux[upper.tri(aux)]
    i <- i+1
  }
  return(x)
}

## Auxiliary for Cerqueira method. The test itself.
c.test <- function(g, gp)
{
  wstat <- sum(abs(colMeans(g)-colMeans(gp)))
  return(wstat)
}

# Auxiliary for Fraiman method. Fraiman test itself according to article.
f.test <- function(g){
  # we need to create a function to calculate a [it's complicated - Appendix 1.3], so far we are using this
  a <- 1

  #how many sets we have
  m <- length(g)

  #size of each set
  l <- unlist(lapply(g, length))

  # make g upper triangular
  g <- lapply(g, f.upper)

  # matrix of mean matrices Mi's
  M <- f.calcM(g)

  # create a list G with all the graphs
  G <- list()
  for(i in 1:length(g)) G <- append(G, g[[i]])

  # calculates $\bar{d}_G(\mathcal{M}_i)$
  # the distance from each $\mathcal{M}_i$ to the entire set of graphs
  sumDG <- rep(0, length(M))
  for(i in 1:length(M)){
    for(j in 1:length(G)){
      sumDG[i] <- sumDG[i] + sum(abs(G[[j]]-M[[i]]))
    }
    sumDG[i] <- sumDG[i]/length(G)
  }

  # calculates $\bar{d}_{G^i}(\mathcal{M}_i)$
  # the distance from each $\mathcal{M}_i$ to the set i of graphs
  sumDGi <- rep(0, length(M))

  for(i in 1:length(M)){
    for(j in 1:length(g[[i]])){
      sumDGi[i] <- sumDGi[i] + sum(abs(g[[i]][[j]] - M[[i]]))
    }
    sumDGi[i] <- sumDGi[i]/length(g[[i]])
  }

  # calculates the final value of the test (equation 2.3)
  # T := \frac{\sqrt(m)}{a} \sum\limits_{i=1}^m \sqrt(n_i) \left( \frac{n_i}{n_i-1} \bar{d}_{G^i}(\mathcal{M_i}) - \frac{n}{n-1}\bar{d}_G(\mathcal{M}_i) \right)
  t1 <- (l/(l-1))*sumDGi
  t2 <- (sum(l)/(sum(l)-1))*sumDG
  t <- (sqrt(m)/a)*sum(sqrt(l)*(t1-t2))

  return(t)
}

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
f.upper <- function(x) lapply(x, function(s){
  s2 <- s
  s2[lower.tri(s2)]<-0
  eval.parent(substitute(s<-s2))
})

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
f.add <- function(x){ list(Reduce("+", x), length(x)) }

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
f.div <- function(x) { x[[1]]/x[[2]] }

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
f.calcM <- function(x){ mapply(f.div, mapply(f.add, x, SIMPLIFY = F), SIMPLIFY = F)}

## Auxiliary for Fraiman method. Padronize input.
f.transform <- function(g)
{
  if(methods::is(g,"igraph")){ return(as.matrix(igraph::get.adjacency(g))) }
  else if(methods::is(g,"list") && methods::is(g[[1]],"igraph")){
    d <- lapply(g, f.transform)
    return(d)
  }
  else if(methods::is(g,"list") && methods::is(g[[1]],"list")){
    d <- lapply(g, f.transform)
    return(d)
  }
}

## Auxiliary for Fraiman method. Boostrap for the test.
f.sampling.distribution <- function(g, maxBoot = 300)
{

  # creates a list with all the graphs
  G <- list()
  n <- length(g)
  for(i in 1:n) G <- append(G, g[[i]])
  m <- length(G)

  dist.boot = c()
  # bootstrap
  for (i_per in 1:maxBoot){
    G1 <- sample(G, m, replace=F)
    #modification made on April 1s, 2019:
    if(n==2){
      l <- list(G1[1:floor(m/2)], G1[(floor(m/2)+1):m])
    }
    else{
      l <- list(G1[1:floor(m/3)], G1[(floor(m/3)+1):(2*floor(m/3))], G1[((2*floor(m/3))+1):m])
    }
    ## original was:
    # l <- list(G1[1:floor(m/2)], G1[(floor(m/2)+1):m])

    dist.boot[i_per] <- f.test(l)
  }
  return(dist.boot)
}

# Graph models  ----------------------------------------------------------------

# Erdos-Renyi graph
ER <- function(n, p, as_matrix = TRUE) {
  if(as_matrix == TRUE){
    return(as.matrix(igraph::get.adjacency(igraph::sample_gnp(n, p))))
  } else {
    return(igraph::sample_gnp(n, p))
  }
}

# Geometric graph
GRG <- function(n, r, as_matrix = TRUE) {
  if(as_matrix == TRUE){
    return (as.matrix(igraph::get.adjacency(igraph::sample_grg(n, r))))
  } else {
    return (igraph::sample_grg(n, r))
  }
}

# Barabasi-Albert graph
BA <- function(n, ps, M = 1, as_matrix = TRUE) {
  if(as_matrix == TRUE){
    return (as.matrix(igraph::get.adjacency(igraph::sample_pa(n, power = ps, m = M,
                                                              directed = FALSE))))
  } else {
    return (igraph::sample_pa(n, power = ps, m = M, directed = FALSE))
  }
}

# Watts-Strogatz graph
WS <- function(n, pr, K = 8, as_matrix = TRUE) {
  if(as_matrix == TRUE){
    return (as.matrix(igraph::get.adjacency(igraph::sample_smallworld(1, n, K, pr))))
  } else {
    return (igraph::sample_smallworld(1, n, K, pr))
  }
}

# K-regular graph
KR <- function(n, k) {
  return(as.matrix(igraph::get.adjacency(igraph::sample_k_regular(n, k))))
}


# Watts-Strogatz small-world graph
WSfun <- function(K){
  f <- function(n, pr, as_matrix = TRUE) {
    WS(n, pr, K=K, as_matrix = as_matrix)
  }
  return(f)
}

# Barabasi-Albert scale-free graph
BAfun <- function(M){
  f <- function(n, ps, as_matrix = TRUE) {
    BA(n, ps, M=M, as_matrix = as_matrix)
  }
  return(f)
}



# new cost function method
# @eigenvalues -> the normalized eigenvalues of the graph of which we want its parameter
# @f           -> the p.d.f of the model
new_cost <- function(eigenvalues, f){
  sum = 0
  for(lambda in eigenvalues){
    found = -1
    dist = Inf
    for(k in 1:length(f$x)){
      if(abs(lambda - f$x[k]) < dist){
        dist = abs(lambda - f$x[k])
        found = k
      }
    }
    if(found != -1)
      sum = sum + log(f$y[found] + 1e-10)
    else
      sum = sum + log(1e-10)
  }
  return (sum)
}

# GIC with repetitions
k_GIC <- function(A, model, p=NULL, bandwidth="Silverman", eigenvalues=NULL,
                  k = 10) {
  gic_estimator = c()
  for(i  in 1:10){
    gic_estimator = append(gic_estimator,GIC(A,model,p,bandwidth,eigenvalues)$value)
  }
  return (mean(gic_estimator))
}

# L2 distance
distance <- function(f1,f2){
  y <- abs(f1$y - f2$y)
  return (trapezoidSum(f1$x,y))
}

## Auxiliary for Tang method.
t.test.stat <- function(X, Y, sigma) {
  n <- nrow(X)
  m <- nrow(Y)
  tmpXX <- sum(exp(-(as.matrix(stats::dist(X))^2)/(2*sigma^2)))
  tmpYY <- sum(exp(-(as.matrix(stats::dist(Y))^2)/(2*sigma^2)))
  tmpXY <- sum(exp(-(t.rect.dist(X,Y))/(2*sigma^2)))
  tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
  return((m+n)*tmp)
}

## Auxiliary for Tang method.
t.embed.graph <- function(g, dim) {
  defaults = igraph::arpack_defaults
  defaults$maxiter = .Machine$integer.max
  lpv = igraph::embed_adjacency_matrix(g,dim, options = defaults)$X

  # Fix signs of eigenvectors issue
  for (i in 1:dim) {
    if (sign(lpv[1, i]) != 1) {
      lpv[, i] = -lpv[, i]
    }
  }
  return(lpv)
}

## Auxiliary for Tang method.
t.rect.dist <- function(X,Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))

  D <- tmp2 - 2*tmp1 + tmp3
  return(D)
}

## Auxiliary for Tang method.
t.get.sigma <- function(X1, X2) {
  v1 = as.vector(stats::dist(X1))
  v2 = as.vector(stats::dist(X2))
  v = base::append(v1, v2)
  sigma = stats::median(v)
  return(sigma)
}

## Auxiliary for Tang method.
t.sampling.distribution <- function(G1, dim, bootstrap_sample_size) {
  Xhat1 = t.embed.graph(G1,dim)
  P = t(Xhat1)
  test_distribution = c()
  i = 1
  while (i <= bootstrap_sample_size) {
    tryCatch({
      G_a = suppressWarnings(igraph::sample_dot_product(P))
      G_b = suppressWarnings(igraph::sample_dot_product(P))
      Xhat_a = suppressWarnings(t.embed.graph(G_a, dim))
      Xhat_b = suppressWarnings(t.embed.graph(G_b, dim))
      sigma = t.get.sigma(Xhat_a, Xhat_b)
      ts = t.test.stat(Xhat_a, Xhat_b, sigma)
      test_distribution[i] = ts
      i = i + 1
    }, error=function(e) {stop(print(e))})
  }
  test_distribution
}

## Auxiliary for Tang method.
t.p_value <- function(ts, test_distribution) {
  area = sum(test_distribution >= ts) / length(test_distribution)
  return(area)
}

## Auxiliary for Tang method.
t.validateInput <- function(G1, G2, dim, maxBoot) {

  !methods::is(G2,"igraph")
  if (methods::is(G1,"dgCMatrix")) { G1 = igraph::graph_from_adjacency_matrix(G1) }
  if (methods::is(G1,"matrix")) { G1 = igraph::graph_from_adjacency_matrix(G1) }
  if (!methods::is(G1,"igraph")) { stop("Input object 'G1' is not an igraph object.") }
  if (methods::is(G2,"dgCMatrix")) { G2 = igraph::graph_from_adjacency_matrix(G2) }
  if (methods::is(G2,"matrix")) { G2 = igraph::graph_from_adjacency_matrix(G2) }
  if (!methods::is(G2,"igraph")) { stop("Input object 'G2' is not an igraph object.") }
  if (!is.null(dim)) {
    if (!methods::is(dim,"numeric") && !is.integer(dim)) { stop("Input 'dim' is not a number.") }
    if (dim%%1 != 0) { stop("Input 'dim' must be an integer.") }
    if (length(dim) > 1) { stop("Input 'dim' has length > 1.") }
    if (dim < 1) { stop("Number of dimensions 'dim' is less than 1.") }
    if (dim >= igraph::gorder(G1) || dim >= igraph::gorder(G2)) { stop("Num. Embedded dimensions 'dim' is greater or equal than number of vertices.") }
  }


  #if (!methods::is(alpha,"numeric")) {
  #  stop("Input object 'alpha' is not a numeric value.")
  #} else if (length(alpha) != 1) {
  #  stop("Input object 'alpha' is not a numeric value.")
  #} else {
  #  if (alpha >= 1 || alpha <= 0) {
  #    stop("Significance level alpha must be strictly between 0 and 1.")
  #  }
  #}

  if (!methods::is(maxBoot,"numeric")) {
    stop("Input object 'maxBoot' is not a numeric value.")
  } else if (length(maxBoot) != 1) {
    stop("Input object 'maxBoot' is not a numeric value.")
  } else {
    if (maxBoot <= 20) {
      stop("The size of bootstrap sample is too small. Pick a larger value.")
    }
  }
  #if (!is.logical(printResult)) { stop("Error: Input 'printResult' must be a logical.")}
}

# Auxiliary for Ghoshdastidar method. Ghoshdastidar test itself according to article.
g.test <- function(x, y){
  m <- min(length(x),length(y))

  x <- lapply(x, function(s){
    s2 <- s
    s2[lower.tri(s2)]<-0
    eval.parent(substitute(s<-s2))
  })

  y <- lapply(y, function(s){
    s2 <- s
    s2[lower.tri(s2)]<-0
    eval.parent(substitute(s<-s2))
  })

  n <- dim(x[[1]])[1]
  Sm1 <- matrix(0,1,m)
  Sp1 <- matrix(0,1,m)
  Sm2 <- matrix(0,1,m)
  Sp2 <- matrix(0,1,m)

  for(i in 1:m){
    Sm1[i] <- sum(x[[i]]*(x[[i]]%*%x[[i]]))/(6*choose(n,3))
    Sp1[i] <- Sm1[i]*log(n)/choose(n,3)

    Sm2[i] <- sum(y[[i]]*(y[[i]]%*%y[[i]]))/(6*choose(n,3))
    Sp2[i] <- Sm2[i]*log(n)/choose(n,3)
  }

  num <- abs(sum(Sm1-Sm2));
  den <- 2*(sqrt(sum(Sp1))+sqrt(sum(Sp2)))

  if(den == 0) den <- 1
  stat <- abs(num/den)
  return(stat)
}


## Auxiliary for Ghoshdastidar method. Boostrap for the test.
g.sampling.distribution <- function(x, y, maxBoot = 300)
{
  m1 <- length(x)
  m2 <- length(y)
  test_distribution = c()
  for (i_per in 1:maxBoot){
    xe <- sample(append(x,y), m1, replace = TRUE)
    ye <- sample(append(x,y), m2, replace = TRUE)
    test_distribution[i_per] <- g.test(xe, ye)
  }
  return(sort(test_distribution))
}

## Auxiliary for Ghoshdastidar method. Padronize input.
g.transform <- function(g)
{
  if(methods::is(g,"igraph")){
    return(list(as.matrix(igraph::get.adjacency(g))))
  }
  else{
    result <- lapply(g, function(x){
      return(as.matrix(igraph::get.adjacency(x)))
    }
    )}
}


## FROM THIS POINT TAIANE CODE
#=============================

#=========================================
## Clustering functions

#' Clustering Expectation-Maximization for Graphs (graph.cem)
#'
#' \code{graph.cem} clusters graphs following an expectation-maximization algorithm based
#' on the Kullback-Leibler divergence between the spectral densities of the
#' graph and of the random graph model.
#'
#' @param g a list containing the graphs or their adjacency matrices to be
#' clustered.
#'
#' @param  model a string that indicates one of the following random graph
#' models: "ER" (Erdos-Renyi random graph), "GRG" (geometric random graph), "KR"
#' (k regular graph), "WS" (Watts-Strogatz model), and "BA" (Barabasi-Albert
#' model).
#'
#' @param k an integer specifying the number of clusters.
#'
#' @param max_iter the maximum number of expectation-maximization steps to execute.
#'
#' @param ncores the number of cores to be used for the parallel processing. The
#' default value is 1.
#'
#' @param bandwidth string showing which criterion is used to choose the
#' bandwidth during the spectral density estimation. Choose between the
#' following criteria: "Silverman" (default), "Sturges", "bcv", "ucv" and "SJ".
#' "bcv" is an abbreviation of biased cross-validation, while "ucv" means
#' unbiased cross-validation. "SJ"  implements the methods of Sheather & Jones
#' (1991) to select the bandwidth using pilot estimation of derivatives.
#'
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{cluster}{a vector of the same length of \code{g} containing the clusterization
#' labels.}
#' \item{parameters}{a vector containing the estimated parameters for the groups.
#' It has the length equals to \code{k}.}
#'
#' @keywords graph.cem
#'
#' @references
#' Celeux, Gilles, and Gerard Govaert. "Gaussian parsimonious clustering
#' models." Pattern recognition 28.5 (1995): 781-793.
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' _Journal of the Royal Statistical Society series B_, 53, 683-690.
#' http://www.jstor.org/stable/2345597.
#'
#' @examples
#'  set.seed(42)
#'  g <- list()
#'  for(i in 1:2){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=0.5)
#'  }
#'  for(i in 3:4){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=1)
#'  }
#'  res <- graph.cem(g, model="ER", k=2, max_iter=1, ncores=1)
#'  res
#' @export
graph.cem <- function(g, model, k, max_iter = 10, ncores=1,
                      bandwidth="Sturges"){

  data.name <- deparse(substitute(g))

  if(methods::is(g,"list") && methods::is(g[[1]],"igraph")){
    g <- f.transform(g)
  }
  `%dopar%` <- foreach::`%dopar%`
  `%:%` <- foreach::`%:%`

  #doMC::registerDoMC(ncores)
  # parallel
  cl <- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)

  tipo <- model
  tau <- matrix(0, nrow = k, ncol = length(g))
  kl <- matrix(0, nrow = k, ncol = length(g))

  prevlik <- 0
  lik <- 1
  count <- 0
  prevlabels <- array(0, length(g))
  vertices <- nrow(g[[1]])
  labels <- array(0, length(g))
  g_GIC <- array(0, length(g))
  p <- array(0, k)

  ## Range that the parameters will be estimated
  if (tipo == "ER") {
    parameters <- seq(0.1, 1, 0.01)
  }
  if (tipo == "GRG") {
    parameters <- seq(0.1, sqrt(2), 0.01)
  }
  if (tipo == "WS") {
    parameters <- seq(0.01, 1, 0.01)
  }
  if (tipo == "KR") {
    parameters <- as.integer(seq(2, 10, 1))
  }
  if (tipo == "BA") {
    parameters <- seq(0.01, 4, 0.01)
  }

  ## Pre-processing of the graph spectra
  eigenvalues <- list()
  j = 1 # added by Grover
  eigenvalues <- foreach::foreach(j = 1:length(g)) %dopar% {
    as.numeric(eigen(g[[j]], only.values = TRUE,
                     symmetric=TRUE)$values)/sqrt(vertices)
  }

  p_graph <- array(0, length(g))
  # Added by Grover
  list_functions = c("GIC","matchFunction","ER","KR","WS","BA","GRG","WSfun","BAfun","modelSpectralDensity","nDensities","gaussianDensity","kernelBandwidth","trapezoidSum","KL","distance")
  ## Parameter estimation
  ret <- foreach::foreach(i = 1:length(g),.export = c("graph.param.estimator",list_functions)) %dopar% {
    graph.param.estimator(g[[i]], model=tipo, parameters=parameters, bandwidth=bandwidth, eigenvalues=eigenvalues[[i]], eps=0.01)
  }
  #
  for(i in 1:length(g)){
    p_graph[i] <- ret[[i]]$p
    g_GIC[i] <- ret[[i]]$KLD
  }

  #Initialize cluster parameters
  p_uniq <- unique(p_graph)
  for(i in 1:k){
    p[i] <- quantile(p_uniq, i/(k+1))
    #the KR parameter needs to be even
    if(tipo == "KR") p[i] <- round(p[i])
  }

  convergiu <- 0
  count <- 0
  while(!convergiu){
    kl <- foreach::foreach(i = 1:k, .combine=cbind,.export = list_functions) %:% foreach::foreach(j = 1:length(g), .combine=c,.export = list_functions) %dopar% {
      GIC(g[[j]], tipo, p[i], bandwidth = bandwidth, eigenvalues = eigenvalues[[j]], dist="KL" )$value
    }
    kl <- t(kl)
    kl[which(kl == Inf)] <- max(kl[which(kl < Inf)])
    kl[which(kl == 0)] <- 0.000000001

    for(i in 1:length(g)){
      tau[,i] <- (1/kl[,i])/sum(1/kl[,i])
    }

    for(i in 1:length(g)){
      labels[i] <- which(tau[,i] == max(tau[,i]))
    }
    #Check if there is an empty group
    for(i in 1:k){
      if(length(which(labels == i))==0) labels[which(tau[i,] == max(tau[i,]))] <- i
    }

    # Estimates the value of p for the models to maximize tae
    for(i in 1:k){
      p[i] <- sum(p_graph[which(labels==i)])/length(which(labels==i))
      if(tipo == "KR") p[i] <- round(p[i])
    }

    prevlik <- lik
    lik <- sum(tau*kl)
    count <- count + 1
    if(count > max_iter){
      convergiu = 1
    }

    if((prevlik!=0 && prevlik/lik > 0.99 && prevlik/lik < 1.01) ) convergiu <- 1
    prevlabels <- labels
  }
  ret <- list("cluster"=labels, "parameters"=p)
  # close cluster
  parallel::stopCluster(cl)


  ###################################################
  method    <- "Clustering Expectation-Maximization for Graphs "

  info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  output     <- list(method=method, info=info,
                     data.name=data.name,
                     cluster=ret$cluster, parameters=ret$parameters)
  attr(output, "class") <- "statGraph"

  return(output)




}

#====================================
#Kmeans

#' K-means for Graphs
#'
#' \code{graph.kmeans} clusters graphs following a k-means algorithm based on the
#' Jensen-Shannon divergence between the spectral densities of the graphs.
#'
#' @param x a list containing the graphs or their adjacency matrices to be
#' clustered.
#'
#' @param k an integer specifying the number of clusters.
#'
#' @param nstart the number of trials of k-means clusterizations. The algorithm
#' returns the clusterization with the best silhouette.
#'
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{cluster}{a vector of the same length of \code{x} containing the
#' clusterization labels.}
#'
#' @keywords k-means
#'
#' @references
#' MacQueen, James. "Some methods for classification and analysis of
#' multivariate observations." Proceedings of the fifth Berkeley symposium on
#' mathematical statistics and probability. Vol. 1. No. 14. 1967.
#'
#' Lloyd, Stuart. "Least squares quantization in PCM." IEEE transactions on
#' information theory 28.2 (1982): 129-137.
#'
#' @examples
#' set.seed(42)
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
graph.kmeans <- function(x, k, nstart=2) {

  data.name <- deparse(substitute(x))
  if(methods::is(x,"list") && methods::is(x[[1]],"igraph")){
    x <- f.transform(x)
  }
  sil <- -1
  num.graphs <- length(x)
  tmp <- spectral_density(x)
  spectral.density <- tmp$spectral.density

  if (k > nstart) nstart <- k

  for (ns in 1:nstart) {
    ## random initialization of the clusters
    label <- sample(seq(1:k), num.graphs, replace=TRUE)
    converged <- FALSE
    while(converged == FALSE) {
      centroid <- matrix(0, k, 512)
      for (j in 1:k) {
        for (i in 1:512) {
          centroid[j,i] <- mean(spectral.density[which(label==j),i])
        }
        centroid[j,] <- centroid[j,] / trapezoidSum(tmp$x, centroid[j,])
      }

      distance <- matrix(0, num.graphs, k)
      for (j in 1:k) {
        tmp1 <- list()
        tmp1$y <- centroid[j,]
        tmp1$x <- tmp$x
        for(i in 1:num.graphs) {
          tmp2 <- list()
          tmp2$y <- spectral.density[i,]
          tmp2$x <- tmp$x
          distance[i,j] <- sqrt(JS(tmp1, tmp2))
        }
      }

      label.new <- array(0, num.graphs)
      for(i in 1:num.graphs) {
        label.new[i] <- which(distance[i,] == min(distance[i,]))[1]
      }
      i <- 1
      while(i<=k) {
        if(length(which(label.new == i)) != 0) {
          i <- i + 1
        }
        else { ## there is an empty cluster
          size.cluster <- array(0,k)
          for(j in 1:k) {
            size.cluster[j] <- length(which(label.new == j))
          }
          largest.cluster <- which(size.cluster == max(size.cluster))
          item <- which(distance[, largest.cluster] ==
                          max(distance[which(label.new == largest.cluster), largest.cluster]))
          label.new[item] <- i
          i <- 1
        }
      }

      if(length(which(label == label.new)) == num.graphs) {
        converged <- TRUE
        sil.new <- mean(cluster::silhouette(label, distance_matrix(tmp))[,3])

        if(sil.new > sil) {
          sil <- sil.new
          label.final <- label
        }
      }
      label <- label.new
    }
  }

  ###################################################
  method    <- "K-means for Graphs"

  info     <- "Clustering the graphs following a k-means algorithm"
  #info     <- paste("Considering a total of ",k," clusters",sep='')

  value     <- list(method=method, info=info,
                    data.name=data.name, cluster=label.final)
  attr(value, "class") <- "statGraph"

  return(value)
}

#Kmeans auxiliary functions
spectral_density <- function(x) {
  num.graphs <- length(x)
  spectrum <- list()

  spectrum[[1]] <- eigen(x[[1]])$values
  max.value <- max(spectrum[[1]])
  min.value <- min(spectrum[[1]])
  for(i in 2:num.graphs) {
    spectrum[[i]] <- eigen(x[[i]])$values
    if(max(spectrum[[i]]) > max.value) {
      max.value <- max(spectrum[[i]])
    }
    if(min(spectrum[[i]]) < min.value) {
      min.value <- min(spectrum[[i]])
    }
  }

  spectral.density <- matrix(0, num.graphs, 512)
  for (i in 1:num.graphs) {
    bw = SturgesBandwidth(spectrum[[i]])
    tmp <- density(spectrum[[i]], bw=bw, from=min.value, to=max.value)
    area <- trapezoidSum(tmp$x, tmp$y)
    spectral.density[i,] <- tmp$y/area ## normalize to AUC == 1
  }
  res <- list()
  res$spectral.density <- spectral.density
  res$x <- tmp$x
  return(res)

}

distance_matrix <- function(x) {
  num.graphs <- nrow(x$spectral.density)
  distance <- matrix(0, num.graphs, num.graphs)
  for (i in 1:(num.graphs-1)) {
    for (j in (i+1):num.graphs) {
      tmp1 <- list()
      tmp1$y <- x$spectral.density[i,]
      tmp1$x <- x$x
      tmp2 <- list()
      tmp2$y <- x$spectral.density[j,]
      tmp2$x <- x$x
      distance[i,j] <- distance[j,i] <- sqrt(JS(tmp1, tmp2))
    }
  }
  return(distance)
}


SturgesBandwidth <- function(x) {
  n <- length(x)
  # Sturges' criterion
  nbins <- ceiling(log2(n) + 1)
  return(abs(max(x) - min(x))/nbins)
}

# FROM THIS POINT GROVER CODE
#############################
#' Degree-based eigenvalue probability
#'
#' \code{fast.eigenvalue.probability} returns the probability of an eigenvalue
#' given the degree and excess degree probability.
#'
#' @param deg_prob The degree probability of the graph.
#'
#' @param q_prob The excess degree probability of the graph.
#'
#' @param all_k  List of sorted unique degrees greater than 1 of the graph.
#'
#' @param z  Complex number whose real part is the eigenvalue whose probability
#' we want to obtain, the imaginary part is a small value (e.g., 1e-3).
#'
#' @param n_iter The maximum number of iterations to perform.
#'
#' @return A complex number whose imaginary part absolute value corresponds to
#' the probability of the given eigenvalue.
#'
#' @keywords eigenvalue_probability
#'
#' @references
#' Newman, M. E. J., Zhang, X., & Nadakuditi, R. R. (2019).
#' Spectra of random networks with arbitrary degrees.
#' Physical Review E, 99(4), 042309.
#'
#' @examples
#' set.seed(42)
#' G <- igraph::sample_smallworld(dim = 1, size = 10, nei = 2, p = 0.2)
#'
#' # Obtain the degree distribution
#' deg_prob <- c(igraph::degree_distribution(graph = G, mode = "all"),0.0)
#' k_deg <- seq(1,length(deg_prob)) - 1
#'
#' # Obtain the excess degree distribution
#' c <- sum(k_deg * deg_prob)
#' q_prob <- c()
#' for(k in 0:(length(deg_prob) - 1)){
#'   aux_q <- (k + 1) * deg_prob[k + 1]/c
#'   q_prob <- c(q_prob,aux_q)
#' }
#'
#' # Obtain the sorted unique degrees greater than 1
#' all_k <- c(1:length(q_prob))
#' valid_idx <- q_prob != 0
#' q_prob <- q_prob[valid_idx]
#' all_k <- all_k[valid_idx]
#'
#' # Obtain the probability of the eigenvalue 0
#' z <- 0 + 0.01*1i
#' eigenval_prob <- -Im(fast.eigenvalue.probability(deg_prob,q_prob,all_k,z))
#' eigenval_prob
#'
#' @export
fast.eigenvalue.probability <- function(deg_prob,q_prob,all_k,z,n_iter = 5000){
  h_z   <- 0 + 0i
  eps <- 1e-7
  all_k_mo <- all_k - 1
  while(n_iter > 0){
    new_h_z <- sum(q_prob/(1 - all_k_mo*h_z))
    new_h_z <- new_h_z/z^2
    # replaces H_z using the new value found
    if(abs(h_z - new_h_z) < eps){
      h_z <- new_h_z
      break
    }
    h_z <- new_h_z
    n_iter <- n_iter - 1
  }

  # returns the final result
  count_z <- 0
  for(k in 1:length(deg_prob)){
    count_z <- count_z + (deg_prob[k])/(1 - (k - 1)*h_z)
  }

  count_z <- (count_z/z)

  return (count_z)
}

#' Degree-based spectral density
#'
#' \code{fast.spectral.density} returns the degree-based spectral density in
#' the interval <\code{from},\code{to}> by using npoints discretization points.
#'
#' @param G The undirected unweighted graph (igraph type) whose spectral
#' density we want to obtain.
#'
#' @param from Lower end of the interval that contain the eigenvalues or
#' smallest eigenvalue of the adjacency matrix of the graph. The smallest
#' eigenvalue is used if the value is not given.
#'
#' @param to  Upper end of the interval that contain the eigenvalues or largest
#' eigenvalue of the adjacency matrix of the graph. The largest eigenvalue is
#' used if the value is not given.
#'
#' @param npoints Number of discretization points of the interval <\code{from},\code{to}>.
#'
#' @param numCores Number of cores to use for parallelization.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{x}{x-value of the the degree-based spectral density of the graph.}
#' \item{y}{y-value of the the degree-based spectral density of the graph.}
#'
#' @keywords eigenvalue_density
#'
#' @references
#' Newman, M. E. J., Zhang, X., & Nadakuditi, R. R. (2019).
#' Spectra of random networks with arbitrary degrees.
#' Physical Review E, 99(4), 042309.
#'
#' @examples
#' set.seed(42)
#' G <- igraph::sample_smallworld(dim = 1, size = 100, nei = 2, p = 0.2)
#'
#' # Obtain the degree-based spectral density
#' density <- fast.spectral.density(G = G, npoints = 80, numCores = 1)
#' density
#'
#' @export
fast.spectral.density <- function(G, from = NULL, to = NULL, npoints = 2000,
                                  numCores = 1){
  graph <- G
  `%dopar%` <- foreach::`%dopar%`
  # Number of vertices
  n <- igraph::vcount(graph)
  # Adjacency matrix
  A <- NULL
  # If 'from' or 'to' are null, then get the adjacency matrix
  if(is.null(from) || is.null(to)){
    A <- igraph::as_adjacency_matrix(graph,type = "both")
  }
  # Obtain the largest eigenvalue
  if(is.null(to)){
    to   <- rARPACK::eigs_sym(A,k = 1)$values[1]
  }
  # Obtain the smallest eigenvalue
  if(is.null(from)){
    from <- rARPACK::eigs_sym(A,k = 1,which = "SA")$values[1]
  }
  # Discretizise interval <\code{from},\code{to}> in npoints
  bw <- (to - from)/npoints
  x <- seq(from,to,bw)
  y <- rep(0,length(x))
  # Obtain the degree and excess degree distribution
  deg_prob <- c(igraph::degree_distribution(graph = graph, mode = "all"),0.0)#/vcount(graph)
  k_deg <- seq(1,length(deg_prob)) - 1
  c <- sum(k_deg * deg_prob)
  q_prob <- c()

  for(k in 0:(length(deg_prob) - 1)){
    aux_q <- (k + 1) * deg_prob[k + 1]/c
    q_prob <- c(q_prob,aux_q)
  }
  # Obtain sorted unique degrees of the graph
  all_k <- c(1:length(q_prob)) #- 1
  valid_idx <- q_prob != 0
  q_prob <- q_prob[valid_idx]
  all_k <- all_k[valid_idx]

  # Obtain the eigenvalue density for each discretized points by using numCores
  # cores.
  #doMC::registerDoMC(numCores)
  cl <- parallel::makePSOCKcluster(numCores)
  doParallel::registerDoParallel(cl)
  i <- NULL
  y <- foreach::foreach(i=1:length(x),.combine = c,.export = c("fast.eigenvalue.probability")) %dopar% {
    z <- x[i] + 0.01*1i
    -Im(fast.eigenvalue.probability(deg_prob,q_prob,all_k,z))
  }

  # close cluster
  parallel::stopCluster(cl)

  ###################################################
  method    <- "Degree-based spectral density"

  info     <- paste("Using",npoints,"discretization points")

  data.name <- deparse(substitute(G))

  output     <- list(method=method, info=info,
                     x=x, y=y)
  attr(output, "class") <- "statGraph"

  return(output)
}

#' Degree-based graph parameter estimator
#'
#' \code{fast.graph.param.estimator} estimates the parameter of the complex
#' network model using the degree-based spectral density and ternary search.
#'
#' @param G The undirected unweighted graph (igraph type).
#'
#' @param model Either a string or a function:
#'
#' A string that indicates one of the following models: "ER" (Erdos-Renyi random
#' graph model), "GRG" (geometric random graph model), "WS" (Watts-Strogatz
#' model), and "BA" (Barabasi-Albert model).
#'
#' A function that returns a Graph generated by a graph model. It must contain
#' two arguments: the first one corresponds to the graph size and the second to
#' the parameter of the model.
#'
#' @param lo Smallest parameter value that the graph model can take.
#'
#' If ``model'' is an string, then the default value of 0 is used for the
#' predefined models ("ER", "GRG", "WS", and "BA").
#'
#' @param hi Largest parameter value that the graph model can take.
#'
#' If ``model'' is an string, then the default values are used for the
#' predefined models 1 for "ER", sqrt(2) for "GRG", 1 for "WS", and 5 for "BA").
#'
#' @param eps Desired precision of the parameter estimate.
#'
#' @param from Lower end of the interval that contain the eigenvalues to
#' generate the degree-based spectral densities. The smallest eigenvalue of the
#' adjacency matrix corresponding to ``graph'' is used if the value is not given.
#'
#' @param to  Upper end of the interval that contain the eigenvalues to generate
#' the degree-based spectral densities. The largest eigenvalue of the adjacency
#' matrix corresponding to ``graph'' is used if the value is not given.
#'
#' @param npoints Number of points to discretize the interval <\code{from},\code{to}>.
#'
#' @param numCores Number of cores to use for parallelization.
#'
#' @return A list with class "statGraph" containing the following components:
#' \item{method}{a string indicating the used method.}
#' \item{info}{a string showing details about the method.}
#' \item{data.name}{a string with the data's name(s).}
#' \item{param}{The degree-based parameter estimate. For the "ER", "GRG", "WS",
#' and "BA" models, the parameter corresponds to the probability to connect a
#' pair of vertices, the radius used to construct the geometric graph in a unit
#' square, the probability to reconnect a vertex, and the scaling exponent
#' respectively.}
#' \item{L1_dist}{The L1 distance between the observed graph and the graph model
#' with the estimated value.}
#'
#' @keywords degree_based_parameter_estimation
#'
#' @examples
#' set.seed(42)
#'
#' ### Example giving only the name of the model to use
#' G <- igraph::sample_smallworld(dim = 1, size = 15, nei = 2, p = 0.2)
#'
#' # Obtain the parameter of the WS model
#' estimated.parameter1 <- fast.graph.param.estimator(G, "WS", lo = 0.1, hi = 0.5,
#'                                                   eps = 1e-1, npoints = 10,
#'                                                   numCores = 1)
#' estimated.parameter1
#'
#' \dontrun{
#' ### Example giving a function instead of a model
#'
#' # Defining the model to use
#' G <- igraph::sample_smallworld(dim = 1, size = 5000, nei = 2, p = 0.2)
#' K <- as.integer(igraph::ecount(G)/igraph::vcount(G))
#' fun_WS <- function(n, param, nei = K){
#'  return (igraph::sample_smallworld(dim = 1,size = n, nei = nei, p = param))
#' }
#'
#' # Obtain the parameter of the WS model
#' estimated.parameter2 <- fast.graph.param.estimator(G, fun_WS, lo = 0.0, hi = 1.0,
#'                                                    npoints = 100, numCores = 2)
#' estimated.parameter2
#' }
#'
#' @export
fast.graph.param.estimator <- function(G, model, lo = NULL, hi = NULL,
                                       eps = 1e-3, from = NULL, to = NULL,
                                       npoints = 2000, numCores = 1){
  graph <- G
  # When the model is a function then check if the smallest and largest values
  # that the parameter can take was also provided

  if(methods::is(model,"function") && (is.null(lo) || is.null(hi))){
    stop("You must specify the largest and smallest parameter value that the graph model can take")
  }
  # If 'model' is a character then define the parameter interval search.
  # Also recover the functions that generate each model
  fun <- model


  if(methods::is(model,"character")){
    if(model == "ER"){
      if(is.null(lo))
        lo <- 0
      if(is.null(hi))
        hi <- 1
      fun <- matchFunction(model)
    } else if(model == "GRG"){
      if(is.null(lo))
        lo <- 0
      if(is.null(hi))
        hi <- sqrt(2)
      fun <- matchFunction(model)
    } else if(model == "WS"){
      if(is.null(lo))
        lo <- 0
      if(is.null(hi))
        hi <- 1
      fun <- WSfun(as.integer(igraph::ecount(graph)/(igraph::vcount(graph))))
    } else if(model == "BA"){
      if(is.null(lo))
        lo <- 0
      if(is.null(hi))
        hi <- 3
      fun <- BAfun(as.integer(igraph::ecount(graph)/(igraph::vcount(graph))))
    } else {
      stop("The 'model' that you specified is not allowed, the allowed model are: \"ER\",\"GRG\",\"WS\" or \"BA\"")
    }
  }

  # Adjacency matrix
  A <- NULL
  # If 'from' or 'up' are null then get the adjacency matrix
  if(is.null(from) || is.null(to)){
    A <- igraph::as_adjacency_matrix(graph,type = "both")
  }
  # Obtain largest eigenvalue
  if(is.null(to)){
    to <- rARPACK::eigs_sym(A,k = 1)$values[1]
  }
  # Obtain smallest eigenvalue
  if(is.null(from)){
    from <- rARPACK::eigs_sym(A,k = 1,which = "SA")$values[1]
  }
  # Obtain the number of vertices of the graph
  n <- igraph::vcount(graph)
  # Obtain density function of the observed graph
  observed_graph_density <- fast.spectral.density(graph, from = from, to = to,
                                                  npoints = npoints,
                                                  numCores = numCores)
  # Make ternary search
  dist_to_1 <- NULL
  dist_to_2 <- NULL
  while(abs(lo - hi) > eps && lo < hi){
    mid1 <- (2*lo + hi)/3
    mid2 <- (2*hi + lo)/3
    # Generate graph using parameter mid1
    if(methods::is(model,"function")){
      Graph1 <- fun(n,mid1)
    } else {
      Graph1 <- fun(n,mid1, as_matrix = FALSE)
    }
    density1 <- fast.spectral.density(Graph1, from = from, to = to,
                                      npoints = npoints,
                                      numCores = numCores)
    # Free memory allocated by Graph1
    rm(Graph1)
    # Generate graph using parameter mid2
    if(methods::is(model,"function")){
      Graph2 <- fun(n,mid2)
    } else {
      Graph2 <- fun(n,mid2, as_matrix = FALSE)
    }
    density2 <- fast.spectral.density(Graph2, from = from, to = to,
                                      npoints = npoints,
                                      numCores = numCores)
    # Free memory allocated by Graph2
    rm(Graph2)
    # Reduce the interval search
    dist_to_1 <- trapezoidSum(observed_graph_density$x,abs(observed_graph_density$y - density1$y))
    dist_to_2 <- trapezoidSum(observed_graph_density$x,abs(observed_graph_density$y - density2$y))
    if(dist_to_1 < dist_to_2){
      hi <- mid2
    } else {
      lo <- mid1
    }
  }
  # The search was finished, now save the parameter and the distance
  param = (lo + hi)/2
  dist  = (dist_to_1 + dist_to_2)/2

  ########################################
  ###################################################
  method    <- "Degree-based graph parameter estimator"

  info     <- "Estimating the parameter of the complex network model using the degree-based spectral density and ternary search"
  #info     <- paste("Using ",bandwidth,"'s criterion to estimate the bandwidth",sep='')

  data.name <- deparse(substitute(G))

  value     <- list(method=method, info=info,
                    data.name=data.name, param=param, L1_dist=dist)
  attr(value, "class") <- "statGraph"

  return(value)
}
