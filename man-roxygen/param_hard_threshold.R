#' @param hard_threshold If specified, the following hard thresholding will be applied
#' to the weights: w = w*(abs(w)>t) where t is the threshold_hard parameter.
#' Defaults to -inf, which means that no threshold will be applied.
#' If all expert weights are thresholded to 0, a weight of 1 will be assigned
#' to the expert with the highest weight prior to thresholding. Thus hard_threshold = 1
#' leads to the 'follow the leader' strategy if method is set to "ewa".
