#' Plot method for online models
#'
#' Plots the most recent weights in each quantile.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @importFrom graphics matplot legend grid polygon
#' @importFrom grDevices rainbow col2rgb rgb
#' @importFrom stats ts.plot
#' @export
plot.online <- function(x, ...) {
    dx <- c(
        dim(x$specification$data$experts)[1],
        dim(x$specification$data$experts[[1]])
    ) # TxDxPxK
    enames <- dimnames(x$experts_loss)[[4]]

    # Univariate point forecasts
    if (dx[2] == 1 && dx[3] == 1) {
        w <- abind::adrop(x$weights, c(2, 3))
        w <- t(apply(w, 1, cumsum))
        w <- cbind(0, w)
        idx <- seq_len(nrow(w))
        cols <- darken(rainbow(n = dim(w)[2]), 1.2)
        plot(idx, idx * NA,
            ylim = c(0, 1),
            ylab = "Weights",
            xlab = "",
            main = "Weights of the Experts",
            xlim = range(idx)
        )
        grid()
        for (i in 2:dim(w)[2]) {
            polygon(
                c(idx, rev(idx)),
                c(w[, i - 1], rev(w[, i])),
                col = cols[i],
                border = cols[i],
                lwd = 1
            )
        }
        legend(
            0, 1,
            legend = enames,
            col = cols[-1],
            lwd = 2,
            cex = 1,
        )
    } else if (dx[2] == 1 && dx[3] != 1) {
        # Univariate probabilistic forecasts
        w <- abind::adrop(x$weights, 2)[nrow(x$weights), , ]
        w <- t(apply(w, 1, cumsum))
        w <- cbind(0, w)
        idx <- seq_len(nrow(w)) / (nrow(w) + 1)
        cols <- darken(rainbow(n = dim(w)[2]), 1.2)
        plot(idx, idx * NA,
            ylim = c(0, 1),
            ylab = "Weights",
            xlab = "Probability",
            main = "Weights of the Experts",
            xlim = range(idx)
        )
        grid()
        for (i in 2:dim(w)[2]) {
            polygon(
                c(idx, rev(idx)),
                c(w[, i - 1], rev(w[, i])),
                col = cols[i],
                border = cols[i],
                lwd = 1
            )
        }
        legend(
            min(idx), 1,
            legend = enames,
            col = cols[-1],
            lwd = 2,
            cex = 1,
        )
    }

    readline(prompt = "Showing Plot 1/2. Press [enter] to continue")

    loss <- cbind(
        "Combination" = apply(x$forecaster_loss, 1, mean),
        apply(x$experts_loss, c(1, 4), mean)
    )
    loss <- apply(loss, 2, cumsum)
    ts.plot(loss,
        xlab = "",
        main = "Cumulative Loss",
        col = cols,
        lwd = 2
    )
    grid()
    legend(
        0, max(loss),
        legend = colnames(loss),
        col = cols, lwd = 2, cex = 1
    )
}