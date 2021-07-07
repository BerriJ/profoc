#' Create a complete ggplot appropriate to a particular data type
#'
#' `autoplot()` uses ggplot2 to draw a particular plot for an object of a
#' particular class in a single command. This defines the S3 generic that
#' other classes and packages can extend.
#'
#' @param object an object, whose class will determine the behaviour of autoplot
#' @param ... other arguments passed to specific methods
#' @return a ggplot object
#' @seealso [autolayer()], [ggplot()] and [fortify()]
#' @export
autoplot <- function(object, ...) {
    UseMethod("autoplot")
}

print_common <- function(x) {
    k <- dim(x$weights)[3]
    experts_loss <- experts <- paste0(
        "E", 1:k, ": ",
        round(apply(x$experts_loss, 3, mean), 4)
    )
    forecasters_loss <- round(mean(x$forecaster_loss), 4)
    cat(
        "Losses: \n",
        experts, "\n",
        "FC:",
        forecasters_loss
    )
}

plot_common <- function(x) {
    weights <- x$weights[nrow(x$weights), , ]
    k <- ncol(weights)
    matplot(
        y = weights,
        x = x$specification$data$tau,
        type = "l",
        ylab = "w",
        xlab = "p",
        ylim = c(0, 1),
        main = "Most Recent Combination Weights",
        lty = 1,
        lwd = 2,
        col = rainbow(k, v = 0.85)
    )

    legend("top", paste("Expert", 1:k),
        bty = "n",
        lty = 1,
        col = rainbow(k, v = 0.85),
        lwd = 2
    )
}

autoplot_common <- function(object) {
    if ("ggplot2" %in% installed.packages()) {
        weights <- object$weights[nrow(object$weights), , ]
        p <- object$specification$data$tau
        weight <- matrix(weights)
        Expert <- as.character(rep(seq_len(ncol(weights)),
            each = nrow(weights)
        ))
        df <- data.frame(weight, Expert, p)
        ggplot2::ggplot(df, ggplot2::aes(x = p, y = weight, fill = Expert)) +
            ggplot2::theme_minimal() +
            ggplot2::geom_area() +
            ggplot2::ggtitle("Most Recent Combination Weights")
    } else {
        cat("Package ggplo2 needs to be installed to use autoplot.")
    }
}