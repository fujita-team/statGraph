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

set.seed(1)
G <- igraph::sample_gnp(n=50, p=0.5)

# Using a string to indicate the graph model
result1 <- GIC(G, "ER", 0.5)
result1

# Using a function to describe the graph model
# Erdos-Renyi graph
model <- function(n, p) {
   return (igraph::sample_gnp(n, p))
}
result2 <- GIC(G, model, 0.5)
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


set.seed(1)
g1 <- g2 <- g3 <- list()
for (i in 1:20) {
  g1[[i]] <- igraph::sample_gnp(50, 0.50)
  g2[[i]] <- igraph::sample_gnp(50, 0.50)
  g3[[i]] <- igraph::sample_gnp(50, 0.52)
}
G <- c(g1, g2, g3)
label <- c(rep(1,20),rep(2,20),rep(3,20))
result <- anogva(G, label, maxBoot=50)
result




cleanEx()
nameEx("graph.acf")
### * graph.acf

flush(stderr()); flush(stdout())

### Name: graph.acf
### Title: Auto Correlation Function Estimation for Graphs
### Aliases: graph.acf
### Keywords: autocorrelation

### ** Examples

set.seed(1)
G <- list()
p <- array(0, 100)
p[1:3] <- rnorm(3)
for (t in 4:100) {
  p[t] <- 0.5*p[t-3] + rnorm(1)
}
ma <- max(p)
mi <- min(p)
p <- (p - mi)/(ma-mi)
for (t in 1:100) {
  G[[t]] <- igraph::sample_gnp(100, p[t])
}
graph.acf(G, plot=TRUE)




cleanEx()
nameEx("graph.cem")
### * graph.cem

flush(stderr()); flush(stdout())

### Name: graph.cem
### Title: Clustering Expectation-Maximization for Graphs (graph.cem)
### Aliases: graph.cem
### Keywords: graph.cem

### ** Examples




cleanEx()
nameEx("graph.cerqueira.test")
### * graph.cerqueira.test

flush(stderr()); flush(stdout())

### Name: graph.cerqueira.test
### Title: Cerqueira's hypothesis testing for random graphs.
### Aliases: graph.cerqueira.test

### ** Examples





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

set.seed(1)
G1 <- G2 <- list()

p <- MASS::mvrnorm(50, mu=c(0,0), Sigma=matrix(c(1, 0.5, 0.5, 1), 2, 2))

ma <- max(p)
mi <- min(p)
p[,1] <- (p[,1] - mi)/(ma - mi)
p[,2] <- (p[,2] - mi)/(ma - mi)

for (i in 1:50) {
  G1[[i]] <- igraph::sample_gnp(50, p[i,1])
  G2[[i]] <- igraph::sample_gnp(50, p[i,2])
}
graph.cor.test(G1, G2)




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
entropy <- graph.entropy(Graph = G)
entropy




cleanEx()
nameEx("graph.fraiman.test")
### * graph.fraiman.test

flush(stderr()); flush(stdout())

### Name: graph.fraiman.test
### Title: Fraiman's hypothesis testing for random graphs.
### Aliases: graph.fraiman.test

### ** Examples





cleanEx()
nameEx("graph.ghoshdastidar.test")
### * graph.ghoshdastidar.test

flush(stderr()); flush(stdout())

### Name: graph.ghoshdastidar.test
### Title: Ghoshdastidar hypothesis testing for large random graphs.
### Aliases: graph.ghoshdastidar.test

### ** Examples





cleanEx()
nameEx("graph.hclust")
### * graph.hclust

flush(stderr()); flush(stdout())

### Name: graph.hclust
### Title: Hierarchical cluster analysis on a list of graphs.
### Aliases: graph.hclust
### Keywords: clustering

### ** Examples

set.seed(1)
G <- list()
for (i in 1:5) {
  G[[i]] <- igraph::sample_gnp(50, 0.5)
}
for (i in 6:10) {
  G[[i]] <- igraph::sample_smallworld(1, 50, 8, 0.2)
}
for (i in 11:15) {
  G[[i]] <- igraph::sample_pa(50, power = 1, directed = FALSE)
}
graph.hclust(G, 3)




cleanEx()
nameEx("graph.kmeans")
### * graph.kmeans

flush(stderr()); flush(stdout())

### Name: graph.kmeans
### Title: K-means for Graphs
### Aliases: graph.kmeans
### Keywords: k-means

### ** Examples

set.seed(42)
g <- list()
for(i in 1:5){
  g[[i]] <- igraph::sample_gnp(30, p=0.2)
}
for(i in 6:10){
  g[[i]] <- igraph::sample_gnp(30, p=0.5)
}
res <- graph.kmeans(g, k=2, nstart=2)
res




cleanEx()
nameEx("graph.model.selection")
### * graph.model.selection

flush(stderr()); flush(stdout())

### Name: graph.model.selection
### Title: Graph model selection
### Aliases: graph.model.selection
### Keywords: model_selection

### ** Examples


## Example using an igraph object as input data
set.seed(1)
G <- igraph::sample_gnp(n=30, p=0.5)

# Using strings to indicate the graph models
result1 <- graph.model.selection(G, models=c("ER", "WS"), eps = 0.5)
result1





cleanEx()
nameEx("graph.mult.scaling")
### * graph.mult.scaling

flush(stderr()); flush(stdout())

### Name: graph.mult.scaling
### Title: Multidimensional scaling of graphs
### Aliases: graph.mult.scaling
### Keywords: multidimensional_scaling

### ** Examples

set.seed(1)
G <- list()
for (i in 1:5) {
  G[[i]] <- igraph::sample_gnp(50, 0.5)
}
for (i in 6:10) {
  G[[i]] <- igraph::sample_smallworld(1, 50, 8, 0.2)
}
for (i in 11:15) {
  G[[i]] <- igraph::sample_pa(50, power = 1, directed = FALSE)
}
graph.mult.scaling(G)




cleanEx()
nameEx("graph.param.estimator")
### * graph.param.estimator

flush(stderr()); flush(stdout())

### Name: graph.param.estimator
### Title: Graph parameter estimator
### Aliases: graph.param.estimator
### Keywords: parameter_estimation

### ** Examples

set.seed(1)
G <- igraph::sample_gnp(n=50, p=0.5)

# Using a string to indicate the graph model
result1 <- graph.param.estimator(G, "ER", eps=0.25)
result1






cleanEx()
nameEx("graph.spectral.density")
### * graph.spectral.density

flush(stderr()); flush(stdout())

### Name: graph.spectral.density
### Title: Graph spectral density
### Aliases: graph.spectral.density
### Keywords: eigenvalue_density

### ** Examples

set.seed(42)
G <- igraph::sample_smallworld(dim = 1, size = 50, nei = 2, p = 0.2)

# Obtain the spectral density
density <- graph.spectral.density(Graph = G)
density




cleanEx()
nameEx("graph.takahashi.test")
### * graph.takahashi.test

flush(stderr()); flush(stdout())

### Name: graph.takahashi.test
### Title: Test for the Jensen-Shannon divergence between graphs
### Aliases: graph.takahashi.test
### Keywords: graph_comparison

### ** Examples

set.seed(1)
G1 <- G2 <- list()
for (i in 1:20) {
  G1[[i]] <- igraph::sample_gnp(n=50, p=0.5)
}
for (i in 1:20) {
  G2[[i]] <- igraph::sample_gnp(n=50, p=0.51)
}
result <- graph.takahashi.test(G1, G2, maxBoot=100)
result




cleanEx()
nameEx("graph.tang.test")
### * graph.tang.test

flush(stderr()); flush(stdout())

### Name: graph.tang.test
### Title: Tang hypothesis testing for random graphs.
### Aliases: graph.tang.test

### ** Examples

set.seed(42)

## test under H0
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
Graph1 <- igraph::sample_dot_product(lpvs)
Graph2 <- igraph::sample_dot_product(lpvs)
D <- graph.tang.test(Graph1,Graph2, 5, printResult = TRUE)

## test under H1
lpvs2 <- matrix(pnorm(200), 20, 10)
lpvs2 <- apply(lpvs2, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g2 <- suppressWarnings(igraph::sample_dot_product(lpvs2))
D <- graph.tang.test(Graph1,Graph2, 5, printResult = TRUE)





cleanEx()
nameEx("sp.anogva")
### * sp.anogva

flush(stderr()); flush(stdout())

### Name: sp.anogva
### Title: Semi-Parametric Analysis Of Graph Variability (ANOGVA)
### Aliases: sp.anogva
### Keywords: semi_parametric_analysis_of_graph_variability

### ** Examples






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
