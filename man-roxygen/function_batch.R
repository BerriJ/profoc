#' Probabilistic Forecast Combination - Batch
#'
#' Returns predictions and weights calculated by sequential numeric
#' optimization. The optimization is done stepwise, always
#' calculating a one-step-ahead forecast.
#'
#' @details batch selects various parameters automatically based on
#' the past loss. For this, the parameters smoothing parameters (see below)
#' can be specified as numeric vectors containing values to consider.
#'
#' This package offers two options for smoothing (Basis Smoothing
#' and P-Splines).
#' Parameters \code{b_smooth} and \code{p_smooth} take named lists to
#' create the corresponding basis and hat matrices. The arguments are:
#' \code{knots} which determines the number of knots to be created, \code{mu},
#' \code{sigma}, \code{sigma}, \code{nonc}, \code{tailweight} correspond to
#' to parameters of the beta distribution, which defines how the knots are
#' #distributed (see \code{?make_knots} for details) the defaults will create
#' an equidistant knot sequence, \code{deg} sets the degree of the spline
#' function and also influences how many outer knots will be used and
#' \code{periodic} which determines whether the spline basis will be periodic.
#' It's possible to provide vectors of values for each of these parameters.
#' In that case, all parameter combinations will be used to create the
#' respective matrices and all candidates will be considered during
#' online-learning. In addition to the inputs mentioned
#' before  \code{p_smooth} requires \code{ndiff} which determines the degree
#' of differentiation applied to the basis-matrix (can take any value
#' between and including 1 and 2), \code{lambda} which determines the degree
#' of penalization applied to the smoothing, higher values will give
#' smoother weight functions. As for the other parameters, it is possible
#' to provide multiple values.
#'
#' @return Returns weights and corresponding predictions.
#' It is possible to impose a convexity constraint to the
#' weights by setting affine and positive to TRUE.
