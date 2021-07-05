#' @param parametergrid User supplied grid of parameters. If i.e. not
#' all combinations of the input vectors should be used. Must be a matrix
#' with 13 columns (online) or 12 columns batch with the following order:
#' basis_knot_distance, basis_knot_distance_power, basis_deg, forget_regret,
#' soft_threshold, hard_threshold, fixed_share, smooth_lambda,
#' smooth_knot_distance, smooth_knot_distance_power, smooth_deg,
#' smooth_ndiff, gamma.
