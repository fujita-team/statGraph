# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

RunTest("graph.takahashi.test", {
            acceptable_error <- 0.1
            G1 <- G2 <- list()
            for (i in 1:20){
                G1[[i]] <- igraph::sample_gnp(n=50, p=0.5)
                G2[[i]] <- igraph::sample_gnp(n=50, p=0.6)
            }
            result <- graph.takahashi.test(G1, G2, maxBoot=100)
            test_that("graph.takahashi.test is working.", {
                          expect_lt(result$p.value, acceptable_error)
            })
})

RunTest("graph.takahashi.test", {
            acceptable_error <- 0.01
            G1 <- G2 <- list()
            for (i in 1:20){
                G1[[i]] <- igraph::sample_gnp(n=10, p=0.5)
                G2[[i]] <- igraph::sample_gnp(n=100, p=0.6)
            }
            result <- graph.takahashi.test(G1, G2, maxBoot=1000)
            test_that("graph.takahashi.test is working.", {
                          expect_lt(result$p.value, acceptable_error)
            })
})


# RunTest("takahashi.test", {
#             acceptable_error <- 0.1
#             G1 <- G2 <- list()
#             for (i in 1:20){
#                 G1[[i]] <- igraph::sample_gnp(n=50, p=0.5, directed=TRUE)
#                 G2[[i]] <- igraph::sample_gnp(n=50, p=0.3, directed=TRUE)
#             }
#             result <- takahashi.test(G1, G2, maxBoot=100, directed=TRUE, distance="L2", npoints=200)
#             test_that("takahashi.test", {
#                           expect_lt(result$p.value, acceptable_error)
#             })
# })
