#' @exportS3Method print statGraph
print.statGraph <- function(x, ...) {
    if (any(names(x) == "method"))
        cat("\n       ", x$method, "\n\n")

    if (any(names(x) == "info"))
        cat("-", x$info, "\n\n")

    if (any(names(x) == "data.name"))
        cat("data:", x$data.name, "\n")

    if (any(names(x) == "entropy"))
        cat("entropy =", x$entropy, "\n")
    if (any(names(x) == "value"))
        cat("value =", x$value, "\n")
    if (any(names(x) == "param"))
        cat("param =", x$param, "\n")
    if (any(names(x) == "dist"))
        cat("dist =", x$dist, "\n")
    if (any(names(x) == "model"))
        cat("model:", x$model, "\n")
    if (any(names(x) == "estimates")) {
        cat("\nestimates: \n")
        print(x$estimates)
    }
    if (any(names(x) == "values")) {
        cat("values:\n")
        print(x$values)
    }

    if (any(names(x) == "cluster")) {
        cat("\ncluster:", x$cluster, "\n")
    }
    if (any(names(x) == "parameters")) {
        cat("\nparameters:", x$parameters, "\n")
    }

    if (any(names(x) == "x")) {
        cat("x:\n")
        print(x$x)
    }
    if (any(names(x) == "y")) {
        cat("\n y:\n")
        print(x$y)
    }
    cat("\n")
}
