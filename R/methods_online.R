#' Predict method for online models
#' Calculates predictions based on new expert advice.
#' This does not update weights. If new observations are available
#' use update instead. The latter updates and computes predictions.
#' @param object Object of class inheriting from 'online'
#' @param new_experts new expert predictions
#' @param ...  further arguments are ignored
#' @return \code{predict.online} produces an updated model object.
#' @importFrom stats predict
#' @rdname online
#' @export
predict.online <- function(object, new_experts, ...) {
    predict_online(object, new_experts)
}

#' Update method for online models
#' Continues learning using new observations and new expert advice.
#' @param object Object of class inheriting from 'online'
#' @param new_y new observations
#' @param new_experts new expert predictions. This must be left unspecified
#' if the model already contains the expert predictions corresponding to new_y.
#' @param ...  further arguments are ignored
#' @return \code{update.online} produces an updated model object.
#' @importFrom stats update
#' @rdname online
#' @export
update.online <- function(object, new_y, new_experts = as.numeric(c()), ...) {
    update_online(object, new_y, new_experts)
}

#' Plot method for online models
#' Plots the most recent weights in each quantile.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @importFrom graphics matplot legend
#' @importFrom grDevices rainbow
#' @rdname online
#' @export
plot.online <- function(x, ...) {
    weights <- x$weights[nrow(x$weights), , ]
    k <- ncol(weights)
    matplot(
        y = weights,
        x = x$specification$data$tau,
        type = "l",
        ylab = "w",
        xlab = "p",
        ylim = c(0, 1),
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

#' @export
autoplot <- function(object, ...) {
    UseMethod("autoplot")
}

#' @rdname online
#' @export
autoplot.online <- function(object, ...) {
    if ("ggplot2" %in% installed.packages()) {
        weights <- object$weights[nrow(object$weights), , ]
        p <- object$specification$data$tau
        weight <- matrix(weights)
        Expert <- as.character(rep(seq_len(ncol(weights)), each = nrow(weights)))
        df <- data.frame(weight, Expert, p)
        ggplot2::ggplot(df, ggplot2::aes(x = p, y = weight, fill = Expert)) +
            ggplot2::theme_minimal() +
            ggplot2::geom_area()
    } else {
        cat("Package ggplo2 needs to be installed to use autoplot.")
    }
}

#' Print method for online models
#' Prints the average loss of all experts (E) and the forecast combination (FC).
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @rdname online
#' @export
print.online <- function(x, ...) {
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