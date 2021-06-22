#' @param soft_threshold If specified the following soft threshold will be applied
#' to the weights: w = sgn(w)*max(abs(w)-t,0) where t is the soft_threshold parameter.
#' Defaults to -inf which means that no threshold will be applied.
