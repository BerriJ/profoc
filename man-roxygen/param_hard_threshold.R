#' @param hard_threshold If specified the following hard thresholding will be applied
#' to the weights: w = w*(abs(w)>t) where t is the threshold_hard parameter.
#' Defaults to -inf which means that no threshold will be applied.
#' If all expert weights are thresholded to 0 a weight of 1 will be assigned
#' to the expert with highest weights prior to thresholding. Thus soft_threshold = 1
#' leads to 'follow the leader' if method = "ewa".
