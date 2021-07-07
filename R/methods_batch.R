#' Print method for batch models
#' Prints the average loss of all experts (E) and the forecast combination (FC).
#' @param x Object of class inheriting from 'batch'
#' @param ...  further arguments are ignored
#' @rdname batch
#' @export
print.batch <- function(x, ...) {
    print_common(x)
}

#' Plot method for batch models
#' Plots the most recent weights in each quantile.
#' @param x Object of class inheriting from 'batch'
#' @param ...  further arguments are ignored
#' @importFrom graphics matplot legend
#' @rdname batch
#' @export
plot.batch <- function(x, ...) {
    plot_common(x)
}

#' Autoplot method for batch models
#' Plots the most recent weights in each quantile using ggplot2.
#' @param object Object of class inheriting from 'batch'
#' @importFrom utils installed.packages
#' @rdname batch
#' @export
autoplot.batch <- function(object, ...) {
    autoplot_common(object)
}