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