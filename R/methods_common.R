print_common <- function(x) {
    expert_names <- dimnames(x$specification$data$experts)[[3]]

    experts_loss <- round(apply(x$experts_loss, 3, mean), 5)

    forecaster_loss <- mean(x$forecaster_loss)

    print(data.frame(
        Name = c(expert_names, "Combination"),
        Loss = c(experts_loss, forecaster_loss)
    ), row.names = FALSE, right = FALSE)
}

plot_common <- function(x) {
    weights <- x$weights[nrow(x$weights), , ]
    k <- ncol(weights)

    expert_names <- dimnames(x$specification$data$experts)[[3]]

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

    legend("top", expert_names,
        bty = "n",
        lty = 1,
        col = rainbow(k, v = 0.85),
        lwd = 2
    )
}

autoplot_common <- function(object) {
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
}
