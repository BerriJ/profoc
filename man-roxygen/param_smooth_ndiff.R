#' @param smooth_ndiff Degree of the differencing operator in the smoothing equation.
#' 1.5 (default) leads to shrikage towards a constant. Can also be 2 or any value
#' in between. If a value in between is used, a weighted sum of the first and second
#' differentiation matrix is calculated.