#' @param parametergrids User supplied grids of parameters. Can be used if not
#' all combinations of the input vectors should be considered. Must be a named
#' list of five matrices. The matrices in list must be named as: "general",
#' "b_smooth_pr", "b_smooth_mv", "p_smooth_pr", "p_smooth_mv".
#' The "general" matrix must contain 11 named columns:
#' "forget_regret", "soft_threshold", "hard_threshold", "fixed_share",
#' "basis_pr_idx", "basis_mv_idx", "hat_pr_idx", "hat_mv_idx",
#' "gamma", "loss_share", "regret_share".
#' The matrices determining the basis smoothing (b_smooth_pr, b_smooth_mv) must
#' contain the following named columns:
#' n, mu, sigma, nonc, tailw, deg, periodic.
#' In addition to the columns of the basis smoothing matrices, the matrices
#' determining the penalized smoothing (p_smooth_pr, p_smooth_mv) must contain
#' the following columns:
#' diff, lambda.
#' The *_idx columns in the general matrix determine which row of the
#' corresponding smoothing matrix is used.
