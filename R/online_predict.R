#' Predict method for online models
#'
#' Calculates predictions based on new expert advice.
#' This does not update weights. If new observations are available
#' use update instead. The latter updates and weights and computes predictions.
#' @param object Object of class inheriting from 'online'
#' @param new_experts new expert predictions
#' @param ...  further arguments are ignored
#' @return \code{predict.online} produces an updated model object.
#' @importFrom stats predict
#' @export
predict.online <- function(object, new_experts, update_model = TRUE, ...) {
    edim <- dim(new_experts)
    if (length(edim) == 3) {
        if (ncol(object$specification$data$y) > 1) { # multivariate point
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
        } else if (ncol(object$specification$data$y) == 1) { # univariate probabilistic
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
    object$weights <- array_to_list(object$weights)
    object <- predict_online(object, new_experts, update_model)

    if (update_model) {
        object$weights <- list_to_array(object$weights)
        return(object)
    } else {
        return(object[[1]])
    }
}