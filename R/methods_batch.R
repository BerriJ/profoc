#' Print method for batch models
#'
#' Prints the average loss of all and the forecast combination.
#' @param x Object of class inheriting from 'batch'
#' @param ...  further arguments are ignored
#' @export
print.batch <- function(x, ...) {
    print_common(x)
}

#' Plot method for batch models
#'
#' Plots the most recent weights in each quantile.
#' @param x Object of class inheriting from 'batch'
#' @param ...  further arguments are ignored
#' @importFrom graphics matplot legend
#' @export
plot.batch <- function(x, ...) {
    plot_common(x)
}

#' Autoplot method for batch models
#'
#' Plots the most recent weights in each quantile using ggplot2.
#' @param object Object of class inheriting from 'batch'
#' @param ...  further arguments are ignored
#' @importFrom utils installed.packages
#' @export
autoplot.batch <- function(object, ...) {
    autoplot_common(object)
}