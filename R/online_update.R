#' Update method for online models
#'
#' Continues learning using new observations and new expert advice.
#' @param object Object of class inheriting from 'online'
#' @param new_y new observations
#' @param new_experts new expert predictions. This must be left unspecified
#' @param trace If a progress bar shall be shown. Defaults to FALSE
#' if the model already contains the expert predictions corresponding to new_y.
#' @param ...  further arguments are ignored
#' @return `update.online` produces an updated model object.
#' @importFrom stats update
#' @export
update.online <- function(object,
                          new_y,
                          new_experts = NULL,
                          trace = FALSE, ...) {
    y <- object$specification$data$y
    names <- object$specification$data$names
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
                new_experts <- array(
                    unlist(new_experts),
                    dim = c(edim[1], edim[2], 1, edim[3])
                )
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

    model_instance <- new(conline)
    model_instance$trace <- trace
    model_instance$init_update(
        object,
        new_y,
        new_experts
    )
    model_instance$learn()

    object <- list(
        predictions = model_instance$predictions,
        predictions_got_sorted = model_instance$predictions_got_sorted,
        weights = model_instance$weights,
        forecaster_loss = model_instance$loss_for,
        experts_loss = model_instance$loss_exp,
        past_performance = model_instance$past_performance,
        opt_index = model_instance$opt_index + 1, # Respect one-based indexing
        parametergrid = model_instance$params,
        params_basis_pr = model_instance$params_basis_pr,
        params_basis_mv = model_instance$params_basis_mv,
        params_hat_pr = model_instance$params_hat_pr,
        params_hat_mv = model_instance$params_hat_mv
    )

    object[["specification"]] <-
        list(
            data =
                list(
                    y = model_instance$y,
                    experts = model_instance$experts,
                    tau = model_instance$tau
                ),
            objects =
                list(
                    weights_tmp = model_instance$weights_tmp,
                    predictions_grid = model_instance$predictions_grid,
                    cum_performance = model_instance$cum_performance,
                    hat_pr = model_instance$hat_pr,
                    hat_mv = model_instance$hat_mv,
                    basis_pr = model_instance$basis_pr,
                    basis_mv = model_instance$basis_mv,
                    V = model_instance$V,
                    E = model_instance$E,
                    eta = model_instance$eta,
                    R = model_instance$R,
                    beta = model_instance$beta,
                    beta0field = model_instance$beta0field
                ),
            parameters =
                list(
                    lead_time = model_instance$lead_time,
                    loss_function = model_instance$loss_function,
                    loss_parameter = model_instance$loss_parameter,
                    loss_gradient = model_instance$loss_gradient,
                    method = model_instance$method,
                    forget_past_performance = model_instance$forget_past_performance,
                    allow_quantile_crossing = model_instance$allow_quantile_crossing,
                    save_past_performance = model_instance$save_past_performance,
                    save_predictions_grid = model_instance$save_predictions_grid
                )
        )

    attr(object, "class") <- c("online", "list")

    model_instance$teardown()
    rm(model_instance)
    object <- post_process_model(model = object, names = names)

    return(object)
}
