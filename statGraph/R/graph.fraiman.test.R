#' Fraiman's hypothesis testing for random graphs.
#'
#' Daniel Fraiman and Ricardo Fraiman test for network differences between
#' groups with an analysis of variance test (ANOVA).
#'
#' Given a list of graphs, the test verifies if all the subpopulations have the
#' same mean network,
#' under the alternative that at least one subpopulation has a different mean
#' network.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param maxPer integer indicating the number of bootstrap resamples
#' (default is \code{300}).
#'
#' @param alpha the significance level for the test (default is \code{0.05}).
#'
#' @param printResult logical indicating if the test must print the result
#' (default is \code{FALSE}).
#'
#' @return A list with class "htest" containing the following components:
#' \itemize{
#' \item{\code{statistic:}}{ the T-value of the test.}
#' \item{\code{p.value:}}{ the p-value of the test.}
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' }
#'
#' @references
#' Fraiman, Daniel, and Ricardo Fraiman. "An ANOVA approach for statistical
#' comparisons of brain networks",
#' https://www.nature.com/articles/s41598-018-23152-5
#'
#' @examples
#' \donttest{
#' set.seed(1)
#'
#' ## test under H0
#' a <- b <- list()
#' for(i in 1:10){
#'   a[[i]] <- igraph::erdos.renyi.game(50,0.5)
#'   b[[i]] <- igraph::erdos.renyi.game(50,0.5)
#' }
#' Graphs <- list(a,b)
#' k <- graph.fraiman.test(Graphs, printResult = TRUE)
#'
#' ## test under H1
#' a <- b <- list()
#' for(i in 1:10){
#'   a[[i]] <- igraph::erdos.renyi.game(50,0.5)
#'   b[[i]] <- igraph::erdos.renyi.game(50,0.6)
#' }
#' Graphs <- list(a,b)
#' k <- graph.fraiman.test(Graphs, printResult = TRUE)
#' }
#'
#' @export
graph.fraiman.test <- function(Graphs, maxPer = 300, alpha = 0.05, printResult = FALSE){
  if(!valid.input(Graphs)) stop("The input should be a list of igraph objects!")
  data.name <- deparse(substitute(Graphs))
  # transform and verify input
  Graphs <- graphList.to.adjList(Graphs)

  D <- fraiman.test(Graphs)

  test_distribution <- fraiman.sampling.distribution(Graphs, maxPer)

  ## modification made in april 1st, 2019:
  p_val <- mean(test_distribution <= D)

  #htest method
  statistic <- D
  names(statistic) <- "T"
  method_info <- "Test for network differences between groups with an analysis of variance test (ANOVA)"
  rval <- list(statistic=statistic, p.value=p_val, method=method_info, data.name=data.name)
  class(rval) <- "htest"
  return(rval)
}


# Auxiliary for Fraiman method. Fraiman test itself according to article.
fraiman.test <- function(g){
  # we need to create a function to calculate a [it's complicated - Appendix 1.3], so far we are using this
  a <- 1

  #how many sets we have
  m <- length(g)

  #size of each set
  l <- unlist(lapply(g, length))

  # make g upper triangular
  g <- lapply(g, fraiman.upper)

  # matrix of mean matrices Mi's
  M <- fraiman.calcM(g)

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
fraiman.upper <- function(x) lapply(x, function(s){
  s2 <- s
  s2[lower.tri(s2)]<-0
  eval.parent(substitute(s<-s2))
})

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
fraiman.add <- function(x){ list(Reduce("+", x), length(x)) }

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
fraiman.div <- function(x) { x[[1]]/x[[2]] }

# Auxiliary for Fraiman method. Intends to speed calculations using R builtins.
fraiman.calcM <- function(x){ mapply(fraiman.div, mapply(fraiman.add, x, SIMPLIFY = F), SIMPLIFY = F)}

## Auxiliary for Fraiman method. Boostrap for the test.
fraiman.sampling.distribution <- function(Graph, maxPer = 300)
{

  # creates a list with all the graphs
  G <- list()
  n <- length(Graph)
  for(i in 1:n) G <- append(G, Graph[[i]])
  m <- length(G)

  dist.boot = c()
  # bootstrap
  for (i_per in 1:maxPer){
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

    dist.boot[i_per] <- fraiman.test(l)
  }
  return(dist.boot)
}
