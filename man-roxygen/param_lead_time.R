#' @param lead_time offset for expert forecasts. Defaults to 0, which means that
#' experts forecast t+1 at t. Setting this to h means experts predictions refer
#' to t+1+h at time t. The weight updates delay accordingly.
