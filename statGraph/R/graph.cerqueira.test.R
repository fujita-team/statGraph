#' Cerqueira's hypothesis testing for random graphs.
#'
#' Andressa Cerqueira, Daniel Fraiman, Claudia D. Vargas and Florencia Leonardi
#' non-parametric test of hypotheses to verify if two samples of random graphs
#' were originated from the same probability distribution.
#'
#' Given two identically independently distributed (idd) samples of graphs \code{Graphs1} and
#' \code{Graphs2}, the test verifies if they have the same distribution by calculating the
#' mean distance D from \code{Graphs1} to \code{Graphs2}. The test rejects the null hypothesis if D is
#' greater than the \code{(1-alpha)}-quantile of the distribution of the test under the
#' null hypothesis.
#'
#' @param Graphs1 a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param Graphs2 a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param maxPer integer indicating the number of bootstrap resamples (default
#' is \code{300}).
#'
#' @param alpha the significance level for the test (default is \code{0.05}).
#'
#' @param printResult logical indicating if the test must print the result
#' (default is \code{FALSE}).
#'
#' @return A list with class "htest" containing the following components:
#' \itemize{
#' \item{\code{statistic:}}{ the W-value of the test.}
#' \item{\code{p.value:}}{ the p-value of the test.}
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{data.name:}}{ a string with the data's name(s).}
#' }
#'
#' @references
#' Andressa Cerqueira, Daniel Fraiman, Claudia D. Vargas and Florencia Leonardi.
#' "A test of hypotheses for random graph distributions built from EEG data",
#' https://ieeexplore.ieee.org/document/7862892
#'
#' @examples
#' \donttest{
#' require(igraph)
#' set.seed(1)
#'
#' ## test under H0
#' Graphs1 <- Graphs2 <- list()
#' for(i in 1:10){
#'   Graphs1[[i]] <- igraph::erdos.renyi.game(50,0.5)
#'   Graphs2[[i]] <- igraph::erdos.renyi.game(50,0.5)
#' }
#' k <- graph.cerqueira.test(Graphs1, Graphs2, printResult = TRUE)
#'
#' ## test under H1
#' Graphs1 <- Graphs2 <- list()
#' for(i in 1:10){
#'   Graphs1[[i]] <- igraph::erdos.renyi.game(50,0.5)
#'   Graphs2[[i]] <- igraph::erdos.renyi.game(50,0.6)
#' }
#' k <- graph.cerqueira.test(Graphs1, Graphs2, printResult = TRUE)
#' }
#'
#' @export
graph.cerqueira.test <- function(Graphs1, Graphs2, maxPer = 300, alpha = 0.05, printResult = FALSE)
{
  if(!valid.input(Graphs1) || !valid.input(Graphs2)) stop("Parameters must be an a list of igraph objects!")

  data.name <- paste(deparse(substitute(Graphs1)), "and", deparse(substitute(Graphs2)))

  Graphs1 <- cerqueira.transform(Graphs1)
  Graphs2 <- cerqueira.transform(Graphs2)

  D <- cerqueira.test(Graphs1,Graphs2)

  test_distribution <- cerqueira.sampling.distribution(Graphs1,Graphs2,maxPer)

  p_val <- mean(test_distribution >= D)

  reject_threshold = quantile(test_distribution,1-alpha)

  if (printResult) {
    if( D > reject_threshold){
      msg <- paste("Reject the null hypothesis that two graphs are identically",
                   "distributed.")
      print(msg)
    } else {
      msg <- paste("Fail to reject the null hypothesis that two graphs are",
                   "identically distributed.")
      print(msg)
    }
  }

  #htest method
  statistic <- D
  names(statistic) <- "W"
  method_info <- "Verify if two samples of random graphs were originated from the same probability distribution."
  rval <- list(statistic=statistic, p.value=p_val, method=method_info, data.name=data.name)
  class(rval) <- "htest"
  return(rval)
}


## Auxiliary for Cerqueira method. Test distribution under the null hypothesis
cerqueira.sampling.distribution <- function(g, gp, maxPer = 300)
{

  m <- nrow(g)+nrow(gp)
  test_distribution = c()
  for (i_per in 1:maxPer){
    total <- rbind(g, gp)
    ind <- sample(1:m, floor(m/2), replace=F)
    xa <- total[ind,]
    ya <- total[-ind,]
    test_distribution[i_per] <- cerqueira.test(xa,ya)
  }
  return(sort(test_distribution))
}

## Auxiliary for Cerqueira method. Fix input format.
cerqueira.transform <- function(Graphs, n = igraph::vcount(Graphs[[1]]))
{
  x <- matrix(0, length(Graphs), n*(n-1)/2)
  i <- 1
  for(graph in Graphs){
    aux <- as.matrix(igraph::get.adjacency(graph))
    x[i,] <- aux[upper.tri(aux)]
    i <- i+1
  }
  return(x)
}

## Auxiliary for Cerqueira method. The test itself.
cerqueira.test <- function(x, y)
{
  wstat <- sum(abs(colMeans(x)-colMeans(y)))
  return(wstat)
}
