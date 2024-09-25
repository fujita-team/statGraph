#' @exportS3Method summary statGraph
summary.statGraph <- function(object, ...) {
    x = object
    if (any(names(x) == "method"))
        cat("method: ", x$method, "\n")

    if (any(names(x) == "info"))
        cat("info: ", x$info, "\n")

    if (any(names(x) == "data.name"))
        cat("data.name: ", x$data.name, "\n")
    if (any(names(x) == "n"))
        cat("Number of objects: ", x$n, "\n")
    if (any(names(x) == "entropy"))
        cat("entropy: ", x$entropy, "\n")
    if (any(names(x) == "value"))
        cat("value: ", x$value, "\n")
    if (any(names(x) == "param"))
        cat("param: ", x$param, "\n")
    if (any(names(x) == "dist"))
        cat("dist: ", x$dist, "\n")
    if (any(names(x) == "model"))
        cat("model: ", x$model, "\n")
    if (any(names(x) == "estimates")) {
        cat("estimates:\n")
        print(x$estimates)
    }
    if (any(names(x) == "values")) {
        cat("values:\n")
        print(x$values)
    }
    if (any(names(x) == "cluster")) {
        cat("cluster: ", x$cluster, "\n")
    }
    if (any(names(x) == "parameters")) {
        cat("parameters: ", x$parameters, "\n")
    }
    if (any(names(x) == "x")) {
        cat("x: \n")
        print(x$x)
    }
    if (any(names(x) == "y")) {
        cat("\n y:\n")
        print(x$y)
    }
}
