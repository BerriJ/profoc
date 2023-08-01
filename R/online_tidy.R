#' @importFrom tidyr tibble
#' @importFrom generics tidy
#' @export
tidy.online.weights <- function(x, ...) {
    names <- x$specification$data$names

    dimnames(x$weights) <- list(
        t = 1:dim(results$weights)[1],
        d = names$y,
        p = probs,
        k = names$experts[[4]]
    )

    return(x)
}
