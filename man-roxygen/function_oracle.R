#' Probabilistic Forecast Combination - Oracle
#'
#' Returns predictions and weights calculated by numeric
#' optimization. The optimization is done in hindsight.
#' This means all observations are used.
#'
#' @return Returns weights and corresponding predictions.
#' It is possible to calculate the best convex combination
#' of weights by setting affine and positive to TRUE.
