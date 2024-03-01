#' Hierarchical cluster analysis on a list of graphs.
#'
#' Given a list of graphs, \code{graph.hclust} builds a hierarchy of clusters
#' according to the Jensen-Shannon divergence between graphs.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param k the number of clusters.
#'
#' @param clus_method the agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of '"ward.D"', '"ward.D2"', '"single"',
#' '"complete"', '"average"' (= UPGMA), '"mcquitty"' (= WPGMA), '"median"'
#' (= WPGMC) or '"centroid"' (= UPGMC).
#'
#' @param dist string indicating if you want to use the "JS" (default), "L1" or "L2"
#' distances. "JS" means Jensen-Shannon divergence.
#'
#' @param ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{hclust:}}{ an object of class \code{hclust} which describes the tree produced
#' by the clustering process.}
#' \item{\code{cluster:}}{ the clustering labels for each graph.}
#' }
#'
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
#' @import methods
#' @export
graph.hclust <- function(Graphs, k, clus_method="complete", dist = "JS", ...) {

  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")

  Graphs <- set.list.spectral.density(Graphs,...)

  d <- distance_matrix(Graphs,dist = dist)

  tmp <- hclust(as.dist(d), method=clus_method)

  res <- list()
  res$hclust <- tmp
  res$cluster <- cutree(tmp, k)

  return(res)
}


#====================================
#Kmeans

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
#' @param dist string indicating if you want to use the "JS" (default), "L1" or "L2"
#' distances. "JS" means Jensen-Shannon divergence.
#'
#' @param  ... Other relevant parameters for \code{\link{graph.spectral.density}}.
#'
#' @return a vector of the same length of \code{Graphs} containing the clusterization
#' labels.
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
graph.kmeans <- function(Graphs, k, nstart=2,dist = "JS",...) {
  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")

  nGraphs <- length(Graphs)
  Graphs <- set.list.spectral.density(Graphs,...)

  sil <- -1

  if (k > nstart) nstart <- k

  for (ns in 1:nstart) {
    ## random initialization of the clusters
    label <- sample(seq(1:k), nGraphs, replace=TRUE)

    converged <- FALSE
    while(converged == FALSE) {
      centroid <- list()
      for (j in 1:k) {
        centroid[[j]] <- get.mean.spectral.density(Graphs[label == j])
      }

      distance_mat <- matrix(0, nGraphs, k)
      for (j in 1:k) {
        for(i in 1:nGraphs) {
          distance_mat[i,j] <- distance(Graphs[[i]]$density,centroid[[j]],dist = dist)
        }
      }

      label.new <- array(0, nGraphs)
      for(i in 1:nGraphs) {
        label.new[i] <- which(distance_mat[i,] == min(distance_mat[i,]))[1]
      }
      i <- 1
      while(i <= k) {
        if(length(which(label.new == i)) != 0) {
          i <- i + 1
        }
        else { ## there is an empty cluster
          size.cluster <- array(0,k)
          for(j in 1:k) {
            size.cluster[j] <- length(which(label.new == j))
          }
          largest.cluster <- which(size.cluster == max(size.cluster))
          item <- which(distance_mat[, largest.cluster] ==
                          max(distance_mat[which(label.new == largest.cluster), largest.cluster]))
          label.new[item] <- i
          i <- 1
        }
      }

      if(sum(label == label.new) == nGraphs) {
        converged <- TRUE
        sil.new <- mean(cluster::silhouette(label, distance_matrix(Graphs,dist = dist))[,3])

        if(sil.new > sil) {
          sil <- sil.new
          label.final <- label
        }
      }
      label <- label.new
    }
  }
  return(label.final)
}

#' Clustering Expectation-Maximization for Graphs (graph.cem)
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
#' models: "ER" (Erdos-Renyi random graph), "GRG" (geometric random graph), "KR"
#' (k regular graph), "WS" (Watts-Strogatz model), and "BA" (Barabási-Albert
#' model).
#'
#' @param k an integer specifying the number of clusters.
#'
#' @param max_iter the maximum number of expectation-maximization steps to execute.
#'
#' @param ... Other relevant parameters for \code{\link{graph.param.estimator}}.
#'
#' @return a list containing three fields:
#' \itemize{
#' \item{\code{labels:}}{ a vector of the same length of \code{Graphs} containing the clusterization labels.}
#' \item{\code{parameters:}}{ a vector containing the estimated parameters for the groups. It has the
#' length equals to \code{k}.}
#' }
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
#' \donttest{
#'  set.seed(42)
#'  g <- list()
#'  for(i in 1:2){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=0.5)
#'  }
#'  for(i in 3:4){
#'    g[[i]] <- igraph::sample_gnp(n=10, p=1)
#'  }
#'  res <- graph.cem(g, model="ER", k=2, max_iter=1)
#'  res
#'  }
#' @export
graph.cem <- function(Graphs, model, k, max_iter = 10, ...){
  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")
  ## Pre-processing of the graph spectra
  Graphs <- set.list.spectral.density(Graphs, ...)

  nGraphs <- length(Graphs)
  tipo <- model
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
  ret <- Map(f = function (G) { graph.param.estimator(Graphs[[i]],model = model,...) },Graphs)
  #
  for(i in 1:nGraphs){
    p_graph[i] <- ret[[i]]$param
    g_GIC[i] <- ret[[i]]$dist
  }

  #Initialize cluster parameters
  p_uniq <- unique(p_graph)
  for(i in 1:k){
    p[i] <- quantile(p_uniq, i/(k+1))
    #the KR parameter needs to be even
    if(tipo == "KR") p[i] <- round(p[i])
  }

  converged <- 0
  count <- 0
  while(!converged){
    kl <- matrix(0,nrow = k,ncol = nGraphs)
    for(i in 1:k){
      for(j in 1:nGraphs){
        kl[i,j] = GIC(Graph = Graphs[[j]], model = model, p = p[i], ...)
      }
    }

    kl[which(kl == Inf)] <- max(kl[which(kl < Inf)])
    kl[which(kl == 0)] <- 1e-9

    #for(i in 1:nGraphs){
    #  tau[,i] <- (1/kl[,i])/sum(1/kl[,i])
    #}
    kl <- 1/kl
    colsum_kl <- colSums(kl)
    tau <- kl/colsum_kl


    for(i in 1:nGraphs){
      labels[i] <- which(tau[,i] == max(tau[,i]))[1]
    }
    #Check if there is an empty group
    for(i in 1:k){
      if(length(which(labels == i))==0) labels[which(tau[i,] == max(tau[i,]))] <- i
    }

    # Estimates the value of p for the models to maximize O tae
    for(i in 1:k){
      p[i] <- sum(p_graph[which(labels==i)])/length(which(labels==i))
      if(model == "KR") p[i] <- round(p[i])
    }

    prevlik <- lik
    lik <- sum(tau*kl)
    count <- count + 1
    if(count > max_iter){
      converged = TRUE
    }

    if((prevlik!=0 && prevlik/lik > 0.99 && prevlik/lik < 1.01)){
      converged <- TRUE
    }
    prevlabels <- labels
  }
  ret <- list("cluster"=labels, "parameters" = p)

  return(ret)
}




