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
    object <- model_instance$output()
    model_instance$teardown()
    rm(model_instance)
    object <- post_process_model(model = object, names = names)

    return(object)
}
