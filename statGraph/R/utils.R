# obtain the eigenvalues of the graph, if Graph contains eigenvalues as attribute then
# such values are returned
graph.eigenvalues <- function(Graph){
  if(!is.null(Graph$density)){
    return (NULL)
  }
  if(!is.null(Graph$eigenvalues)){
    return (Graph$eigenvalues)
  } else {
    A <- igraph::get.adjacency(Graph)
    eigenvalues <- as.numeric(eigen(A, only.values=TRUE, symmetric=TRUE)$values)
    eigenvalues <- eigenvalues/sqrt(nrow(A))
    rm(A)
    return (eigenvalues)
  }
}

# checks if the variable is a graph, a list of graphs, or a list of lists of graphs
valid.input <- function(Graph){
  if(methods::is(Graph,"list")) return (Reduce(f = "&",Map(f = valid.input,Graph)))
  if(methods::is(Graph,"igraph")) return (TRUE)
  return (FALSE)
}


# Returns the density function for a sample x at n points in the interval [from, to]
gaussianDensity <- function(x, from=NULL, to=NULL, bandwidth="Silverman",
                            npoints=1024) {
  if (bandwidth == "Sturges"){
    bw <- kernelBandwidth(x)
  }else if (bandwidth == "Silverman"){
    bw <- stats::bw.nrd0(x)
  }else if (bandwidth == "bcv"){
    bw <- suppressWarnings(stats::bw.bcv(x))
  }else if (bandwidth == "ucv"){
    bw <- suppressWarnings(stats::bw.ucv(x))
  }else if (bandwidth == "SJ"){
    bw <- "SJ"
  }else{
    stop("Please, choose a valid bandwidth.")
  }
  if (bw == 0){
    # this case happens when all eigenvalues are equal
    return(NA)
  }
  if (is.null(from) || is.null(to)){
    f <- stats::density(x, bw=bw, n=npoints)
  }else{
    f <- stats::density(x, bw=bw, from=from, to=to, n=npoints)
  }
  f$y <- f$y + 1e-12 # we do not want the area to be zero, so we add a very small number
  area <- trapezoidSum(f$x, f$y)
  return(list("x"=f$x, "y"=f$y/area,"from" = min(f$x),"to" = max(f$x),"bw" = f$bw,"method"="exact"))
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

# distance between two density functions
distance <- function(f1,f2,dist = "KL"){
  if(dist == "KL") return (KL(f1,f2))
  if(dist == "L1") return (L2(f1,f2))
  if(dist == "L2") return (L2(f1,f2))
  if(dist == "JS") return (JS(f1,f2))
  # error if the distance parameter does not exists
  stop(paste0(dist," distance measure is not valid. Use: KL,L1,L2, or JS"))
}

# Return the L1 norm between two densities
L1 <- function(f1,f2){
  y <- abs(f1$y - f2$y)
  return (trapezoidSum(f1$x, y))
}


# Returns the L2 norm between two densities
L2 <- function(f1,f2){
  y <- (f1$y - f2$y)^2
  return (trapezoidSum(f1$x,y))
}

# Returns the Kullback-Leibler divergence between two densities
KL <- function(f1, f2) {
  y <- f1$y
  diff_zero <- (y != 0)
  eq_zero <- (f2$y == 0)

  if(sum(diff_zero & eq_zero) != 0) return (Inf)

  y[diff_zero] <- y[diff_zero]*log(y[diff_zero]/f2$y[diff_zero])

  #y <- f1$y
  #n <- length(y)
  #for (i in 1:n) {
  #  if (y[i] != 0 && f2$y[i] == 0){
  #    return (Inf)
  #  }
  #  if (y[i] != 0)
  #    y[i] <- y[i]*log(y[i]/f2$y[i])
  #}
  return (trapezoidSum(f1$x, y))
}

# Returns the Jensen-Shannon divergence between two densities
JS <- function(f1, f2) {
  fm <- f1
  fm$y <- (f1$y + f2$y)/2
  return((KL(f1, fm) + KL(f2, fm))/2)
}

# functions to obtain the smallest and largest eigenvalues of a graph or a list of graphs
get.smallest.eigenvalue <- function(Graphs){
  if(methods::is(Graphs,"igraph")){
    if(is.null(Graphs$eigenvalues)){
      A <- igraph::as_adjacency_matrix(Graphs,type = "both")
      ev <- rARPACK::eigs_sym(A,k = 1,which = "SA")$values[1]
      rm(A)
      return (ev)
    } else {
      return (Graphs$eigenvalues[igraph::vcount(Graphs)])
    }
  } else if(methods::is(Graphs,"list")) {
    return (Reduce(f = "min",Map(f = get.smallest.eigenvalue,Graphs)))
  }
  stop("Input should be a Graph or a list of graphs.")
}

# functions to obtain the largest and largest eigenvalues of a graph or a list of graphs
get.largest.eigenvalue <- function(Graphs){
  if(methods::is(Graphs,"igraph")){
    if(is.null(Graphs$eigenvalues)){
      A <- igraph::as_adjacency_matrix(Graphs,type = "both")
      ev <- rARPACK::eigs_sym(A,k = 1)$values[1]
      rm(A)
      return (ev)
    } else {
      return (Graphs$eigenvalues[1])
    }
  } else if(methods::is(Graphs,"list")) {
    return (Reduce(f = "max",Map(f = get.largest.eigenvalue,Graphs)))
  }
  stop("Input should be a Graph or a list of graphs.")
}

# transform a list of graphs to their respective adjacency matrices
graphList.to.adjList <- function(Graphs){
  if(methods::is(Graphs,'igraph')){ return(igraph::get.adjacency(Graphs,type = "both")) }
  else if(methods::is(Graphs,'list') && methods::is(Graphs[[1]],'igraph')){
    d <- lapply(Graphs, graphList.to.adjList)
    return(d)
  }
  else if(methods::is(Graphs,'list') && methods::is(Graphs[[1]],'list')){
    d <- lapply(Graphs, graphList.to.adjList)
    return(d)
  }
}
