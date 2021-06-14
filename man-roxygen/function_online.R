#' Probabilistic Forecast Combination - Online
#'
#' Returns predictions and weights calculated by online-learning algorithms
#' using CRPS Learning. By default, the weights are calculated by
#' gradient based bernstein online aggregation (BOAG).
#'
#' @return online can tune various parameters automatically based on
#' the past loss. For this, lambda, forget, fixed_share, gamma, ndiff,
#' deg and knot_distance can be specified as numeric vectors containing
#' parameters to consider. online will automatically try all possible
#' combinations of values provide.
