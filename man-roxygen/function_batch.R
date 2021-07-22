#' Probabilistic Forecast Combination - Batch
#'
#' Returns predictions and weights calculated by sequential numeric
#' optimization. The optimization is done stepwise, always
#' calculating a one-step-ahead forecast.
#'
#' @return Returns weights and corresponding predictions.
#' It is possible to impose a convexity constraint to the
#' weights by setting affine and positive to TRUE.
