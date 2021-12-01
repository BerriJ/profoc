#' Print method for online models
#'
#' Prints the average loss of all experts and the forecast combination.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @export
print.online <- function(x, ...) {
    print_common(x)
}

#' Predict method for online models
#'
#' Calculates predictions based on new expert advice.
#' This does not update weights. If new observations are available
#' use update instead. The latter updates and computes predictions.
#' @param object Object of class inheriting from 'online'
#' @param new_experts new expert predictions
#' @param ...  further arguments are ignored
#' @return \code{predict.online} produces an updated model object.
#' @importFrom stats predict
#' @export
predict.online <- function(object, new_experts, ...) {
    predict_online(object, new_experts)
}

#' Update method for online models
#'
#' Continues learning using new observations and new expert advice.
#' @param object Object of class inheriting from 'online'
#' @param new_y new observations
#' @param new_experts new expert predictions. This must be left unspecified
#' @param trace If a progressbar shall be shown. Defaults to FALSE
#' if the model already contains the expert predictions corresponding to new_y.
#' @param ...  further arguments are ignored
#' @return \code{update.online} produces an updated model object.
#' @importFrom stats update
#' @export
update.online <- function(object, new_y, new_experts = NULL, trace = FALSE, ...) {
    if (is.vector(new_y)) {
        new_y <- matrix(new_y)
    }

    if (is.null(new_experts)) {
        new_experts <- list()
        dim(new_experts) <- c(0, 0)
    } else {
        edim <- dim(new_experts)

        if (length(edim) == 3) {
            if (ncol(new_y) > 1) { # multivariate point
                new_experts <- array(unlist(new_experts), dim = c(edim[1], edim[2], 1, edim[3]))
                new_experts <- lapply(seq_len(edim[1]),
                    asub,
                    x = new_experts,
                    dims = 1,
                    drop = TRUE
                )
                dim(new_experts) <- c(edim[1], 1)
            } else if (ncol(new_y) == 1) { # univariate probabilistic
                new_experts <- lapply(seq_len(edim[1]),
                    asub,
                    x = new_experts,
                    dims = 1,
                    drop = FALSE
                )
                dim(new_experts) <- c(edim[1], 1)
            }
        } else if (length(edim) == 4) { # multivariate probabilistic
            new_experts <- lapply(seq_len(edim[1]),
                asub,
                x = new_experts,
                dims = 1,
                drop = TRUE
            )
            dim(new_experts) <- c(edim[1], 1)
        }
    }

    object$weights <- array_to_list(object$weights)
    object$past_performance <- array_to_list(object$past_performance)
    object$experts_loss <- array_to_list(object$experts_loss)

    update_online(object, new_y, new_experts, trace)

    dimnames(object$specification$data$y) <- dimnames(new_y)

    object$weights <- list_to_array(object$weights)
    object$past_performance <- list_to_array(object$past_performance)
    object$experts_loss <- list_to_array(object$experts_loss)

    parnames <- c(
        "basis_knot_distance",
        "basis_knot_distance_power",
        "basis_deg",
        "forget_regret",
        "threshold_soft",
        "threshold_hard",
        "fixed_share",
        "p_smooth_lambda",
        "p_smooth_knot_distance",
        "p_smooth_knot_distance_power",
        "p_smooth_deg",
        "smooth_diff",
        "gamma",
        "loss_share",
        "regret_share",
        "mv_basis_knot_distance",
        "mv_basis_knot_distance_power",
        "mv_basis_deg",
        "mv_p_smooth_lambda",
        "mv_p_smooth_knot_distance",
        "mv_p_smooth_knot_distance_power",
        "mv_p_smooth_deg",
        "mv_p_smooth_ndiff"
    )

    colnames(object$chosen_parameters) <- parnames
    colnames(object$parametergrid) <- parnames

    return(object)
}

#' Plot method for online models
#'
#' Plots the most recent weights in each quantile.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @importFrom graphics matplot legend
#' @importFrom grDevices rainbow
#' @export
plot.online <- function(x, ...) {
    plot_common(x, ...)
}

#' Autoplot method for online models
#'
#' Plots the most recent weights in each quantile using ggplot2.
#' @param object Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @importFrom utils installed.packages
#' @export
autoplot.online <- function(object, ...) {
    autoplot_common(object)
}