# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf


RunTest("anogva", {
            acceptable_error <- 0.1
            g1 <- g2 <- g3 <- list()
            for (i in 1:20){
                g1[[i]] <- igraph::sample_gnp(50, 0.50)
                g2[[i]] <- igraph::sample_gnp(50, 0.50)
                g3[[i]] <- igraph::sample_gnp(50, 0.52)
            }
            G <- c(g1, g2, g3)
            label <- c(rep(1, 20), rep(2, 20), rep(3, 20))
            result <- anogva(G, label, maxBoot=100, directed=FALSE)
            test_that("Anogva with undirected graphs.", {
                expect_lt(result$p.value, acceptable_error)
                })
})

RunTest("anogva", {
            acceptable_error <- 0.1
            g1 <- g2 <- g3 <- list()
            for (i in 1:20){
                g1[[i]] <- igraph::sample_gnp(50, 0.50, directed=TRUE)
                g2[[i]] <- igraph::sample_gnp(50, 0.50, directed=TRUE)
                g3[[i]] <- igraph::sample_gnp(50, 0.32, directed=TRUE)
            }
            G <- c(g1, g2, g3)
            label <- c(rep(1, 20), rep(2, 20), rep(3, 20))
            result <- anogva(G, label, maxBoot=50, directed=TRUE, npoints=200)
            test_that("Anogva with directed graphs.", {
                expect_lt(result$p.value, acceptable_error)
               })
})

# TODO: Should we also test H1? We would need to verify if the distribution of p-values obtained when all groups are the same form a uniform distribution.
# This is probably a little bit too intense for automatic testing.
# For now, we are going to be skipping it!
