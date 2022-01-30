#' Summary method for online models
#'
#' Calculates parameters chosen during optimization and aggregates losses.
#' @param object Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @export
summary.online <- function(object, ...) {
    experts_loss <- round(apply(object$experts_loss, 4, mean), 5)
    forecaster_loss <- mean(object$forecaster_loss)

    pargrid <- object$parametergrid[object$opt_index, ]
    chosen_basis_pr <- object$params_basis_pr[pargrid[, "basis_pr_idx"], ]
    chosen_basis_mv <- object$params_basis_mv[pargrid[, "basis_mv_idx"], ]
    chosen_hat_pr <- object$params_hat_pr[pargrid[, "hat_pr_idx"], ]
    chosen_hat_mv <- object$params_hat_mv[pargrid[, "hat_mv_idx"], ]

    lt <- object$specification$parameters$lead_time
    # Set unused values to NA
    if (lt > 0) {
        pargrid[1:lt, ] <- NA
        chosen_basis_pr[1:lt, ] <- NA
        chosen_basis_mv[1:lt, ] <- NA
        chosen_hat_pr[1:lt, ] <- NA
        chosen_hat_mv[1:lt, ] <- NA
    }

    value <- as.numeric(tail(pargrid, 1))

    parnames <- colnames(object$parametergrid)
    cat("General\n")
    df <- data.frame("Value" = value, row.names = parnames)
    print(df)

    cat("\n")

    cat("Probabilistic Basis Smoothing\n")
    print(chosen_basis_pr[nrow(chosen_basis_pr), ])

    cat("\n")
    cat("Multivariate Basis Smoothing\n")
    print(chosen_basis_mv[nrow(chosen_basis_mv), ])

    cat("\n")

    cat("Probabilistic P-Spline Smoothing\n")
    print(chosen_hat_pr[nrow(chosen_hat_pr), ])

    cat("\n")

    cat("Multivariate P-Spline Smoothing\n")
    print(chosen_hat_mv[nrow(chosen_hat_mv), ])

    return(invisible(list(
        experts_loss = experts_loss,
        forecaster_loss = forecaster_loss,
        chosen_parameters = pargrid,
        chosen_basis_pr = chosen_basis_pr,
        chosen_basis_mv = chosen_basis_mv,
        chosen_hat_pr = chosen_hat_pr,
        chosen_hat_mv = chosen_hat_mv
    )))
}