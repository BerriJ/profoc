#' Predict method for online models
#'
#' Calculates predictions based on new expert advice.
#' This does not update weights. If new observations are available
#' use update instead. The latter updates and weights and computes predictions.
#' @param object Object of class inheriting from 'online'
#' @param new_experts new expert predictions
#' @param update_model Defines whether the model object should be updated or not.
#' If TRUE, new forecaster and expert predictions are appended onto
#' the respective object items. Defaults to TRUE.
#' @param ...  further arguments are ignored
#' @return \code{predict.online} produces an updated model object.
#' @importFrom stats predict
#' @importFrom abind abind
#' @importFrom utils tail
#' @export
predict.online <- function(object, new_experts, update_model = TRUE, ...) {
    edim <- dim(new_experts)
    if (length(edim) == 3) {
        if (ncol(object$specification$data$y) > 1) {
            new_experts <- array(
                unlist(new_experts),
                dim = c(edim[1], edim[2], 1, edim[3])
            )
        } else if (ncol(object$specification$data$y) == 1) {
            new_experts <- array(
                unlist(new_experts),
                dim = c(edim[1], 1, edim[2], edim[3])
            )
        }
    }

    # Use most recent weights
    w <- abind::adrop(tail(object$weights, 1), 1)
    new_predictions <- array(NA, dim = dim(new_experts)[-4])
    for (i in seq_len(nrow(new_experts))) {
        # Predict
        new_predictions[i, , ] <-
            apply(w * abind::adrop(
                new_experts[i, , , , drop = FALSE],
                1
            ), 1:2, sum)

        # Sort
        if (!object$specification$parameters$allow_quantile_crossing) {
            new_predictions[i, , ] <- apply(new_predictions[i, , , drop = FALSE], MARGIN = 2, FUN = sort)
        }
    }

    if (update_model) {
        object$predictions <- abind::abind(
            object$predictions,
            new_predictions,
            along = 1
        )
        object$specification$data$experts <- c(
            object$specification$data$experts,
            array_to_list(new_experts)
        )
        dim(object$specification$data$experts) <- c(
            length(object$specification$data$experts), 1, 1
        )
        return(object)
    } else {
        return(new_predictions)
    }
}