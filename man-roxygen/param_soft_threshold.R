#' @param soft_threshold If specified, the following soft threshold will be applied
#' to the weights: w = sgn(w)*max(abs(w)-t,0) where t is the soft_threshold parameter.
#' Defaults to -inf, which means that no threshold will be applied.
#' If all expert weights are thresholded to 0, a weight of 1 will be assigned
#' to the expert with the highest weights prior to thresholding. Thus soft_threshold = 1
#' leads to the 'follow the leader' strategy if method is set to "ewa".
