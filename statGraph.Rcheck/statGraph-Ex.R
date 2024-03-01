pkgname <- "statGraph"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('statGraph')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GIC")
### * GIC

flush(stderr()); flush(stdout())

### Name: GIC
### Title: Graph Information Criterion (GIC)
### Aliases: GIC
### Keywords: graph_information_criterion

### ** Examples

A <- as.matrix(igraph::get.adjacency(igraph::sample_gnp(n=50, p=0.5)))
# Using a string to indicate the graph model
result1 <- GIC(A, "ER", 0.5)
result1

# Using a function to describe the graph model
# Erdos-Renyi graph
model <- function(n, p) {
   return(as.matrix(igraph::get.adjacency(igraph::sample_gnp(n, p))))
}
result2 <- GIC(A, model, 0.5)
result2



cleanEx()
nameEx("anogva")
### * anogva

flush(stderr()); flush(stdout())

### Name: anogva
### Title: ANOGVA Analysis Of Graph Variability
### Aliases: anogva
### Keywords: analysis_of_graph_variability

### ** Examples


require(igraph)
g1 <- g2 <- g3 <- list()
for (i in 1:20) {
   G1 <- erdos.renyi.game(50, 0.50)
   g1[[i]] <- get.adjacency(G1)
   G2 <- erdos.renyi.game(50, 0.50)
   g2[[i]] <- get.adjacency(G2)
   G3 <- erdos.renyi.game(50, 0.52)
   g3[[i]] <- get.adjacency(G3)
}
g <- c(g1, g2, g3)
label <- c(rep(1,20),rep(2,20),rep(3,20))
result <- anogva(g, label, numBoot=50)
result




cleanEx()
nameEx("cerqueira")
### * cerqueira

flush(stderr()); flush(stdout())

### Name: cerqueira
### Title: Andressa Cerqueira, Daniel Fraiman, Claudia D. Vargas and
###   Florencia Leonardi non-parametric test of hypotheses to verify if two
###   samples of random graphs were originated from the same probability
###   distribution.
### Aliases: cerqueira

### ** Examples

## Not run: 
##D require(igraph)
##D set.seed(42)
##D 
##D ## test under H0
##D a <- b <- list()
##D for(i in 1:10){
##D   a[[i]] <- erdos.renyi.game(50,0.5)
##D   b[[i]] <- erdos.renyi.game(50,0.5)
##D }
##D k <- cerqueira(a, b, printResult = TRUE)
##D 
##D ## test under H1
##D a <- b <- list()
##D for(i in 1:10){
##D   a[[i]] <- erdos.renyi.game(50,0.5)
##D   b[[i]] <- erdos.renyi.game(50,0.6)
##D }
##D k <- cerqueira(a, b, printResult = TRUE)
## End(Not run)




cleanEx()
nameEx("fast.eigenvalue.probability")
### * fast.eigenvalue.probability

flush(stderr()); flush(stdout())

### Name: fast.eigenvalue.probability
### Title: Degree-based eigenvalue probability
### Aliases: fast.eigenvalue.probability
### Keywords: eigenvalue_probability

### ** Examples

G <- igraph::sample_smallworld(dim = 1, size = 10, nei = 2, p = 0.2)
# Obtain the degree distribution
deg_prob <- c(igraph::degree_distribution(graph = G, mode = "all"),0.0)
k_deg <- seq(1,length(deg_prob)) - 1
# Obtain the excess degree distribution
c <- sum(k_deg * deg_prob)
q_prob <- c()
for(k in 0:(length(deg_prob) - 1)){
  aux_q <- (k + 1) * deg_prob[k + 1]/c
  q_prob <- c(q_prob,aux_q)
}
# Obtain the sorted unique degrees greater than 1
all_k <- c(1:length(q_prob))
valid_idx <- q_prob != 0
q_prob <- q_prob[valid_idx]
all_k <- all_k[valid_idx]
# Obtain the probability of the eigenvalue 0
z <- 0 + 0.01*1i
eigenval_prob <- -Im(fast.eigenvalue.probability(deg_prob,q_prob,all_k,z))
eigenval_prob




cleanEx()
nameEx("fast.graph.param.estimator")
### * fast.graph.param.estimator

flush(stderr()); flush(stdout())

### Name: fast.graph.param.estimator
### Title: Degree-based graph parameter estimator
### Aliases: fast.graph.param.estimator
### Keywords: degree_based_parameter_estimation

### ** Examples

### Example giving only the name of the model to use
G <- igraph::sample_smallworld(dim = 1, size = 15, nei = 2, p = 0.2)
# Obtain the parameter of the WS model
estimated.parameter <- fast.graph.param.estimator(G, "WS",lo = 0,hi = 0.5,eps = 1e-1, npoints = 10,
                                                  numCores = 1)
estimated.parameter
### Example giving a function instead of a model (uncomment to execute)
# Defining the model to use
#G <- igraph::sample_smallworld(dim = 1, size = 5000, nei = 2, p = 0.2)
#K <- as.integer(igraph::ecount(G)/igraph::vcount(G))
#fun_WS <- function(n, param, nei = K){
# return (igraph::sample_smallworld(dim = 1,size = n, nei = nei,p = param))
#}
# Obtain the parameter of the WS model
#estimated.parameter <- fast.graph.param.estimator(G, fun_WS, lo = 0.0, hi = 1.0,
#                                                   npoints = 100, numCores = 2)
#estimated.parameter




cleanEx()
nameEx("fast.spectral.density")
### * fast.spectral.density

flush(stderr()); flush(stdout())

### Name: fast.spectral.density
### Title: Degree-based spectral density
### Aliases: fast.spectral.density
### Keywords: eigenvalue_density

### ** Examples

G <- igraph::sample_smallworld(dim = 1, size = 100, nei = 2, p = 0.2)
# Obtain the degree-based spectral density
density <- fast.spectral.density(graph = G, npoints = 80, numCores = 1)
density




cleanEx()
nameEx("fraiman")
### * fraiman

flush(stderr()); flush(stdout())

### Name: fraiman
### Title: Daniel Fraiman and Ricardo Fraiman test for network differences
###   between groups with an analysis of variance test (ANOVA).
### Aliases: fraiman

### ** Examples

## Not run: 
##D require(igraph)
##D set.seed(42)
##D 
##D ## test under H0
##D a <- b <- d <- list()
##D for(i in 1:10){
##D   a[[i]] <- erdos.renyi.game(50,0.5)
##D   b[[i]] <- erdos.renyi.game(50,0.5)
##D }
##D d <- list(a,b)
##D k <- fraiman(d, printResult = TRUE)
##D 
##D ## test under H1
##D a <- b <- d <- list()
##D for(i in 1:10){
##D   a[[i]] <- erdos.renyi.game(50,0.5)
##D   b[[i]] <- erdos.renyi.game(50,0.6)
##D }
##D d <- list(a,b)
##D k <- fraiman(d, printResult = TRUE)
## End(Not run)




cleanEx()
nameEx("gCEM")
### * gCEM

flush(stderr()); flush(stdout())

### Name: gCEM
### Title: Clustering Expectation-Maximization for Graphs (gCEM)
### Aliases: gCEM
### Keywords: gCEM

### ** Examples

require(igraph)
 g <- list()
 for(i in 1:2){
   g[[i]] <- igraph::get.adjacency(igraph::sample_gnp(n=10, p=0.5))
 }
 for(i in 3:4){
   g[[i]] <- igraph::get.adjacency(igraph::sample_gnp(n=10, p=1))
 }
 res <- gCEM(g, model="ER", num_clusters=2, max_iter=1, ncores=1)



cleanEx()
nameEx("ghoshdastidar")
### * ghoshdastidar

flush(stderr()); flush(stdout())

### Name: ghoshdastidar
### Title: Ghoshdastidar hypothesis testing for large random graphs.
### Aliases: ghoshdastidar

### ** Examples

## Not run: 
##D require(igraph)
##D set.seed(42)
##D 
##D ## test for sets with more than one graph each under H0
##D x <- y <- list()
##D for(i in 1:10){
##D   x[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50,0.6)))
##D   y[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50,0.6)))
##D }
##D D <- ghoshdastidar(x, y, printResult = TRUE)
##D 
##D ## test for sets with more than one graph each under H1
##D x <- y <- list()
##D for(i in 1:10){
##D   x[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50,0.6)))
##D   y[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50,0.7)))
##D }
##D D <- ghoshdastidar(x, y, printResult = TRUE)
##D 
##D ## test for sets with only one graph each under H0
##D x <- y <- list()
##D x[[1]] <- erdos.renyi.game(300, 0.6)
##D y[[1]] <- erdos.renyi.game(300, 0.6)
##D D <- ghoshdastidar(x, y, two.sample= TRUE, printResult = TRUE)
##D 
##D ## test for sets with only one graph each under H1
##D x <- y <- list()
##D x[[1]] <- erdos.renyi.game(300, 0.6)
##D y[[1]] <- erdos.renyi.game(300, 0.7)
##D D <- ghoshdastidar(x, y, two.sample= TRUE, printResult = TRUE)
## End(Not run)




cleanEx()
nameEx("graph.acf")
### * graph.acf

flush(stderr()); flush(stdout())

### Name: graph.acf
### Title: Auto Correlation Function Estimation for Graphs
### Aliases: graph.acf
### Keywords: autocorrelation

### ** Examples

require(igraph)
x <- list()
p <- array(0, 100)
p[1:3] <- rnorm(3)
for (t in 4:100) {
    p[t] <- 0.5*p[t-3] + rnorm(1)
}
ma <- max(p)
mi <- min(p)
p <- (p - mi)/(ma-mi)
for (t in 1:100) {
    x[[t]] <- get.adjacency(erdos.renyi.game(100, p[t]))
}
graph.acf(x, plot=TRUE)




cleanEx()
nameEx("graph.cluster")
### * graph.cluster

flush(stderr()); flush(stdout())

### Name: graph.cluster
### Title: Hierarchical cluster analysis on a list of graphs.
### Aliases: graph.cluster
### Keywords: clustering

### ** Examples

require(igraph)
g <- list()
for (i in 1:5) {
    g[[i]] <- as.matrix(get.adjacency(
                        erdos.renyi.game(50, 0.5, type="gnp",
                                         directed = FALSE)))
}
for (i in 6:10) {
    g[[i]] <- as.matrix(get.adjacency(
                        watts.strogatz.game(1, 50, 8, 0.2)))
}
for (i in 11:15) {
    g[[i]] <- as.matrix(get.adjacency(
                        barabasi.game(50, power = 1,
                                      directed = FALSE)))
}
graph.cluster(g, 3)




cleanEx()
nameEx("graph.cor.test")
### * graph.cor.test

flush(stderr()); flush(stdout())

### Name: graph.cor.test
### Title: Test for Association / Correlation Between Paired Samples of
###   Graphs
### Aliases: graph.cor.test
### Keywords: correlation_coefficient

### ** Examples

require(igraph)
x <- list()
y <- list()

p <- MASS::mvrnorm(50, mu=c(0,0), Sigma=matrix(c(1, 0.5, 0.5, 1), 2, 2))

ma <- max(p)
mi <- min(p)
p[,1] <- (p[,1] - mi)/(ma - mi)
p[,2] <- (p[,2] - mi)/(ma - mi)

for (i in 1:50) {
    x[[i]] <- get.adjacency(erdos.renyi.game(50, p[i,1]))
    y[[i]] <- get.adjacency(erdos.renyi.game(50, p[i,2]))
}

graph.cor.test(x, y)




cleanEx()
nameEx("graph.entropy")
### * graph.entropy

flush(stderr()); flush(stdout())

### Name: graph.entropy
### Title: Graph spectral entropy
### Aliases: graph.entropy
### Keywords: spectral_entropy

### ** Examples

G <- igraph::sample_gnp(n=100, p=0.5)
A <- as.matrix(igraph::get.adjacency(G))
entropy <- graph.entropy(A)
entropy




cleanEx()
nameEx("graph.model.selection")
### * graph.model.selection

flush(stderr()); flush(stdout())

### Name: graph.model.selection
### Title: Graph model selection
### Aliases: graph.model.selection
### Keywords: model_selection

### ** Examples


require(igraph)
A <- as.matrix(get.adjacency(erdos.renyi.game(30, p=0.5)))
# Using strings to indicate the graph models
result1 <- graph.model.selection(A, models=c("ER", "WS"), eps=0.5)
result1
# Using functions to describe the graph models
# Erdos-Renyi graph
model1 <- function(n, p) {
   return(as.matrix(get.adjacency(erdos.renyi.game(n, p))))
}
# Watts-Strougatz graph
model2 <- function(n, pr, K=8) {
    return(as.matrix(get.adjacency(watts.strogatz.game(1, n, K, pr))))
}
parameters <- list(seq(0, 1, 0.5), seq(0, 1, 0.5))
result2 <- graph.model.selection(A, list(model1, model2), parameters)
result2



cleanEx()
nameEx("graph.mult.scaling")
### * graph.mult.scaling

flush(stderr()); flush(stdout())

### Name: graph.mult.scaling
### Title: Multidimensional scaling of graphs
### Aliases: graph.mult.scaling
### Keywords: multidimensional_scaling

### ** Examples

require(igraph)
g <- list()
for (i in 1:5) {
    g[[i]] <- as.matrix(get.adjacency(
                        erdos.renyi.game(50, 0.5, type="gnp",
                                         directed = FALSE)))
}
for (i in 6:10) {
    g[[i]] <- as.matrix(get.adjacency(
                        watts.strogatz.game(1, 50, 8, 0.2)))
}
for (i in 11:15) {
    g[[i]] <- as.matrix(get.adjacency(
                        barabasi.game(50, power = 1,
                                      directed = FALSE)))
}
graph.mult.scaling(g)




cleanEx()
nameEx("graph.param.estimator")
### * graph.param.estimator

flush(stderr()); flush(stdout())

### Name: graph.param.estimator
### Title: Graph parameter estimator
### Aliases: graph.param.estimator
### Keywords: parameter_estimation

### ** Examples

require(igraph)
A <- as.matrix(get.adjacency(erdos.renyi.game(50, p=0.5)))

# Using a string to indicate the graph model
result1 <- graph.param.estimator(A, "ER", eps=0.25)
result1

## Using a function to describe the graph model
## Erdos-Renyi graph
# model <- function(n, p) {
#    return(as.matrix(get.adjacency(erdos.renyi.game(n, p))))
# }
# result2 <- graph.param.estimator(A, model,  seq(0.2, 0.8, 0.1))
# result2



cleanEx()
nameEx("graph.test")
### * graph.test

flush(stderr()); flush(stdout())

### Name: graph.test
### Title: Test for the Jensen-Shannon divergence between graphs
### Aliases: graph.test
### Keywords: graph_comparison

### ** Examples

library(igraph)
x <- y <- list()
for (i in 1:20)
   x[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50, p=0.5)))
for (i in 1:20)
   y[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50, p=0.51)))

result <- graph.test(x, y, numBoot=100)
result




cleanEx()
nameEx("kmeans.graph")
### * kmeans.graph

flush(stderr()); flush(stdout())

### Name: kmeans.graph
### Title: K-means for Graphs
### Aliases: kmeans.graph
### Keywords: k-means

### ** Examples

require(igraph)
g <- list()
for(i in 1:5){
  g[[i]] <- get.adjacency(sample_gnp(30, p=0.2))
}
for(i in 6:10){
  g[[i]] <- get.adjacency(sample_gnp(30, p=0.5))
}
res <- kmeans.graph(g, k=2, nstart=2)




cleanEx()
nameEx("sp.anogva")
### * sp.anogva

flush(stderr()); flush(stdout())

### Name: sp.anogva
### Title: Semi-Parametric Analysis Of Graph Variability (ANOGVA)
### Aliases: sp.anogva
### Keywords: semi_parametric_analysis_of_graph_variability

### ** Examples


## Please uncomment the following lines to run an example
# require(igraph)
# set.seed(42)
# model <- "ER"
# graph <- list()

## Under H0
# graph[[1]] <- get.adjacency(erdos.renyi.game(50, 0.5))
# graph[[2]] <- get.adjacency(erdos.renyi.game(50, 0.5))
# graph[[3]] <- get.adjacency(erdos.renyi.game(50, 0.5))
# result <- sp.anogva(graph, model, maxBoot = 300)
# result

## Under H1
# graph[[1]] <- get.adjacency(erdos.renyi.game(50, 0.5))
# graph[[2]] <- get.adjacency(erdos.renyi.game(50, 0.55))
# graph[[3]] <- get.adjacency(erdos.renyi.game(50, 0.5))
# result <- sp.anogva(graph, model, maxBoot = 300)
# result




cleanEx()
nameEx("tang")
### * tang

flush(stderr()); flush(stdout())

### Name: tang
### Title: Tang hypothesis testing for random graphs.
### Aliases: tang

### ** Examples

require(igraph)
set.seed(42)

## test under H0
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g1 <- sample_dot_product(lpvs)
g2 <- sample_dot_product(lpvs)
D <- tang(g1,g2, 5, printResult = TRUE)

## test under H1
lpvs2 <- matrix(pnorm(200), 20, 10)
lpvs2 <- apply(lpvs2, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g2 <- suppressWarnings(sample_dot_product(lpvs2))
D <- tang(g1,g2, 5, printResult = TRUE)





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
