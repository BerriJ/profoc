#' @param basis_knot_distance determines the distance of the knots in the probability
#' basis. Defaults to the follwing sequence:
#' c(2^seq(log(1/(length(tau)+1),2)-1, -1, length=5),1).
#' The best knot_distance is selected in each step according
#' to the past performance.
