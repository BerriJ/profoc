#' Predict method for online models
#' Calculates predictions based on new expert advice.
#' This does not update weights. If new observations are available
#' use update instead. The latter updates and computes predictions.
#' @param object Object of class inheriting from 'online'
#' @param new_experts new expert advices
#' @param ...  further arguments are ignored
#' @return \code{predict.online} produces an updated model object.
#' @importFrom stats predict
#' @rdname online
#' @export
predict.online <- function(object, new_experts, ...) {
    predict_online(object, new_experts)
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
        x = x$specification$tau,
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