# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf

# RunTest("GIC", {
#             acceptable_error <- 0.1
#             G <- igraph::sample_gnp(n=50, p=0.5)
#
#             result1 <- GIC(G, "ER", 0.5)
#             test_that("GIC", {
#                 expect_lt(result1$value, acceptable_error)
#             })
# })
# RunTest("GIC", {
#             model <- function(n, p){
#                 return(igraph::sample_gnp(n, p))
#             }
#             result2 <- GIC(G, model, 0.5)
#             test_that("GIC", {
#                 expect_lt(result1$value, acceptable_error)
#             })
# })
