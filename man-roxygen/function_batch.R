#' Probabilistic Forecast Combination - Batch
#'
#' Returns predictions and weights calculated by sequential numeric
#' optimization. The optimization is done stepwise, allways
#' calculating a one-step-ahead forecast.
#'
#' @return Returns weights and corresponding predictions.
#' It is possible to calculate the best convex combination
#' of weights by setting affine and positive to TRUE.
