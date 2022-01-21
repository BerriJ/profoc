#' Probabilistic Forecast Combination - Online
#'
#' Returns predictions and weights calculated by online-learning algorithms
#' using CRPS Learning.
#'
#'
#' @details online selects various parameters automatically based on
#' the past loss. For this, lambda, forget, fixed_share, gamma, and the
#' smoothing parameters (see below) can be specified as numeric vectors
#' containing values to consider.
#'
#' This package offers 2 optios for smoothing (Basis Smoothing
#' and P-Splines). Both options can be used to apply smoothing the weights
#' over dimension D (covariates) or P (quantiles) or both.
#' Parameters \code{b_smooth_pr} and \code{b_smooth_mv} take named lists to create
#' the corresponding basis matrices. The arguments include are: \code{knots}
#' which determines the number of knots to be created (excluding outer
#' knots which are added additionally), \code{beta_a} and \code{beta_b} correspond
#' to the parameters
#' of the beta distribution which defined how the knots are distributed
#' (see \code{?make_knots} for details) both default to 1 which will create an
#' equidistant knot sequence, \code{deg} sets the degree of the spline function
#' and also influences how many outer knots will be used, \code{outer_knots}
#' determines if outer knots should be created, this defaults to TRUE.
#' It's possible to provide vectors of values for each of these parameters.
#' In that case, all parametercombinations will be used to create the
#' respective matrices and all candidates will be considered during
#' online-learning.
#' Parameters \code{p_smooth_pr} and \code{p_smooth_mv} determine the hat-matrix
#' cretion fot P-Spline smoothing. In addition to the inputs mentioned
#' before, they require to provide \code{ndiff} which determines the degree
#' of differentiation applied to the basis-matrix (can take any value
#' between and including 1 and 2), \code{lambda} which determines the degree
#' of penalization applied to the smoothing, higher values will give
#' smoother weight functions. As for the other parameters, it is possible
#' to provide multiple values.
#'
#' @return Returns weights and corresponding predictions.
