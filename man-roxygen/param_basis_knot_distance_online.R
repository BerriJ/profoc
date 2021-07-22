#' @param basis_knot_distance determines the distance of the knots in the probability
#' basis. Defaults to 1 / (dim(experts)[2] + 1), which means that one knot is created
#' for every quantile. Takes a vector with values >0 where values > .5 correspond to constant
#' weights (only one single knot).
