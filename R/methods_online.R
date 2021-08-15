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
#' if the model already contains the expert predictions corresponding to new_y.
#' @param ...  further arguments are ignored
#' @return \code{update.online} produces an updated model object.
#' @importFrom stats update
#' @export
update.online <- function(object, new_y, new_experts = as.numeric(c()), ...) {
    update_online(object, new_y, new_experts)
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