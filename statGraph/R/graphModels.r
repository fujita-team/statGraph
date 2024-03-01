# return the graph model function
get.graph.model <- function(model,mean_deg = 1){
  if(methods::is(model,"function")) return (model)
  else if (methods::is(model,"character")){
    if(model == "BA") return (BAfun(mean_deg))
    if(model == "WS") return (WSfun(as.integer(mean_deg/2)))
    return (match.fun(model))
  }
  stop("Graph model should be either an string ('ER','WS','BA', or 'GRG') or a function")
}

# Erdos-Renyi graph
ER <- function(n, p) {
  return(igraph::sample_gnp(n, p))
}

# Geometric graph
GRG <- function(n, r) {
  return (igraph::sample_grg(n, r))
}

# Barabasi-Albert graph
BA <- function(n, ps, M = 1) {
  return (igraph::sample_pa(n, power = ps, m = M, directed = FALSE))
}

# Watts-Strogatz graph
WS <- function(n, pr, K = 8) {
  return (igraph::sample_smallworld(1, n, K, pr))
}

# K-regular graph
KR <- function(n, k) {
  return (igraph::sample_k_regular(n, k))
}

# Watts-Strogatz small-world graph
WSfun <- function(K){
  f <- function(n, pr) {
    WS(n, pr, K=K)
  }
  return(f)
}

# Barabasi-Albert scale-free graph
BAfun <- function(M){
  f <- function(n, ps) {
    BA(n, ps, M=M)
  }
  return(f)
}


############################################################
# Methods to obtain the degree/exact spectral density of a graph model#
############################################################

# obtain the graph model spectral density
# model: string ER|WS|BA|GRG or a function
# obtain the exact spectral density of a graph model
# ngraphs: number of graphs to generate to obtain the spectral density
graph.model.spectral.density <- function(model,n,p,mean_deg = 1,ngraphs = 50,method = "diag",...){
  if(methods::is(model,"list")){
    if(!is.null(model$x) && !is.null(model$y)){
      model$from <- min(model$x)
      model$to <- max(model$x)
      return (model)
    }
  }
  parameters = list(...)
  parameters$model = model
  parameters$n = n
  parameters$p = p
  parameters$ngraphs = ngraphs
  parameters$from = parameters$from
  parameters$to = parameters$to
  parameters$mean_deg <- mean_deg
  den_fun <- match.fun(paste0("graph.",method,".model.spectral.density"))
  return (do.call(den_fun,parameters))
}

# returns the degree-based spectral density of a graph model
# extra: aditional parameter for the graph model
graph.fast.model.spectral.density <- function(model,n,p,mean_deg = 1,ngraphs = 50,...){

  # recover spectral density parameters
  spectral_density_param = get.fast.density.parameters(...)
  # recover graph model
  graph_gen <- get.graph.model(model,mean_deg)
  # generate a Graph
  Graph <- graph_gen(n,p)
  spectral_density_param$Graph = Graph
  den_fun <- do.call(graph.fast.spectral.density,spectral_density_param)
  rm(Graph)
  return (den_fun)
}


# returns the exact spectral density of a graph model
graph.diag.model.spectral.density <- function(model,n,p,mean_deg = 1,ngraphs = 50,...){
  # recover spectral density parameters
  spectral_density_param = get.diag.density.parameters(...)
  # recover graph model
  graph_gen <- get.graph.model(model,mean_deg)
  # generate ngraphs with graph model
  Graphs <- list()
  for(idx in 1:ngraphs){
    Graph <- graph_gen(n,p)
    Graph$eigenvalues <- graph.eigenvalues(Graph)
    Graphs[[idx]] <- Graph
    rm(Graph)
  }
  # Now obtain the spectral density of each graph
  Graphs <- set.list.spectral.density(Graphs,method = "diag",...)

  return (get.mean.spectral.density(Graphs))
}



# Obtain the smallest eigenvalue
get.model.interval <- function(n,m,model,parameter,eps,search){
  if (methods::is(model,"function")) {
    if(is.null(parameter)){
      stop("It is necessary to enter the parameters that will be evaluated.")
    }
  }
  #
  if(methods::is(parameter,"atomicVector"))
    return (parameter)
  if(methods::is(parameter,"list"))
    return (list(parameter))
  #
  if(search == "grid"){
    parameters <- NULL
    if(is.null(parameter)){
      if (model == "GRG") parameters <- seq(0.05,sqrt(2), eps)
      else if (model == "BA") parameters <- seq(0.05,5, eps)
      else if (model == "KR") parameters <- as.integer(seq(0.05, 1, eps)*n)
      else parameters <- seq(0.05, 1, eps)
      return (parameters)
    } else {
      if(methods::is(parameter,"list")){
        parameters <- seq(parameter$lo,parameter$hi, eps)
      } else if(methods::is(parameter,"atomicVector")) {
        parameters <- parameter
      }
      return (parameters)
    }
  } else {
    intervals <- list(parameter)
    if(methods::is(model,"character")){
      if(model == "ER"){
        p = 2*(m/(n*(n - 1)))
        intervals = list(list("lo" = p,"hi" = p))
      } else if(model == "KR"){
        p = (2*m)/n
        intervals = list(list("lo" = p,"hi" = p))
      } else {
        if (model == "WS") intervals = list(list("lo" = 0.05,"hi" = 0.4),list("lo" = 0.4,"hi" = 1))
        else if (model == "BA") intervals = list(list("lo" = 0.05,"hi" = 1.25),list("lo" = 1.25,"hi" = 2.25),list("lo" = 2.25,"hi" = 3))
        else if (model == "GRG") intervals = list(list("lo" = 0.05,"hi" = 0.8),list("lo" = 0.8,"hi" = 1.4))
      }
    }
    return (intervals)
  }

}
