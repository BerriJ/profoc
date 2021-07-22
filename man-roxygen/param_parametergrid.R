#' @param parametergrid User supplied grid of parameters. Can be used if not
#' all combinations of the input vectors should be considered. Must be a matrix
#' with 13 columns (online) or 12 columns batch with the following order:
#' basis_knot_distance, basis_knot_distance_power, basis_deg, forget_regret,
#' soft_threshold, hard_threshold, fixed_share, p_smooth_lambda,
#' p_smooth_knot_distance, p_smooth_knot_distance_power, p_smooth_deg,
#' p_smooth_ndiff, gamma.
