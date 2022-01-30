#' Autoplot method for online models
#'
#' Plots the most recent weights in each quantile using ggplot2.
#' @param object Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @importFrom utils installed.packages
#' @export
autoplot.online <- function(object, ...) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        weights <- object$weights[nrow(object$weights), , ]
        p <- object$specification$data$tau
        weight <- matrix(weights)
        expert_names <- dimnames(object$specification$data$experts)[[3]]
        Expert <- as.character(rep(expert_names,
            each = nrow(weights)
        ))
        df <- data.frame(weight, Expert, p)
        ggplot2::ggplot(df, ggplot2::aes(x = p, y = weight, fill = Expert)) +
            ggplot2::theme_minimal() +
            ggplot2::geom_area() +
            ggplot2::ggtitle("Most Recent Combination Weights")
    } else {
        cat("Package ggplot2 needs to be installed to use autoplot.")
    }
}