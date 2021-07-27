#' @template function_batch
#'
#' @template param_y
#' @template param_experts
#' @template param_tau
#' @template param_affine
#' @template param_positive
#' @template param_intercept
#' @template param_debias
#' @template param_lead_time
#'
#' @param initial_window Defines the size of the initial estimation window.
#' @param rolling_window Defines the size of the rolling window. Defaults to
#' the value of initial_window. Set it to the number of observations to receive
#' an expanding window.
#'
#' @template param_loss_function
#' @template param_loss_parameter
#'
#' @template param_basis_knot_distance_batch
#' @template param_basis_knot_distance_power
#' @template param_basis_deg_batch
#'
#' @template param_forget
#'
#' @template param_soft_threshold
#' @template param_hard_threshold
#'
#' @template param_fixed_share
#'
#' @template param_p_smooth_lambda
#' @template param_p_smooth_knot_distance
#' @template param_p_smooth_knot_distance_power
#' @template param_p_smooth_deg
#' @template param_p_smooth_ndiff
#'
#' @template param_parametergrid_max_combinations
#' @template param_parametergrid
#' @template param_forget_past_performance
#'
#' @template param_allow_quantile_crossing
#'
#' @template param_trace
#'
#' @examples
#' T <- 50 # Observations
#' N <- 2 # Experts
#' P <- 9 # Quantiles
#' prob_grid <- 1:P / (P + 1)
#'
#' y <- rnorm(n = T) # Realized
#' experts <- array(dim = c(T, P, N)) # Predictions
#' for (t in 1:T) {
#'     experts[t, , 1] <- qnorm(prob_grid, mean = -1, sd = 1)
#'     experts[t, , 2] <- qnorm(prob_grid, mean = 3, sd = sqrt(4))
#' }
#'
#' model <- batch(
#'     y = matrix(y),
#'     experts = experts,
#'     p_smooth_lambda = 10
#' )
#'
#' print(model)
#' plot(model)
#' autoplot(model)
#' @export
batch <- function(y,
                  experts,
                  tau = 1:dim(experts)[2] / (dim(experts)[2] + 1),
                  affine = FALSE,
                  positive = FALSE,
                  intercept = FALSE,
                  debias = TRUE,
                  lead_time = 0,
                  initial_window = 30,
                  rolling_window = initial_window,
                  loss_function = "quantile",
                  loss_parameter = 1,
                  basis_knot_distance = 1 / (dim(experts)[2] + 1),
                  basis_knot_distance_power = 1,
                  basis_deg = 1,
                  forget = 0,
                  soft_threshold = -Inf,
                  hard_threshold = -Inf,
                  fixed_share = 0,
                  p_smooth_lambda = -Inf,
                  p_smooth_knot_distance = basis_knot_distance,
                  p_smooth_knot_distance_power = basis_knot_distance_power,
                  p_smooth_deg = basis_deg,
                  p_smooth_ndiff = 1.5,
                  parametergrid_max_combinations = 100,
                  parametergrid = NULL,
                  forget_past_performance = 0,
                  allow_quantile_crossing = FALSE,
                  trace = TRUE) {

    # Ensure that batch_rcpp does not expand a grid for basis_knot_distance
    # and p_smooth_knot_distance etc.
    if (missing(p_smooth_knot_distance)) {
        p_smooth_knot_distance <- as.numeric(c())
    }

    if (missing(p_smooth_knot_distance_power)) {
        p_smooth_knot_distance_power <- as.numeric(c())
    }

    if (missing(p_smooth_deg)) {
        p_smooth_deg <- as.numeric(c())
    }

    if (is.null(parametergrid)) {
        parametergrid <- matrix(ncol = 0, nrow = 0)
    }

    model <- batch_rcpp(
        y = y,
        experts = experts,
        tau = tau,
        affine = affine,
        positive = positive,
        intercept = intercept,
        debias = debias,
        lead_time = lead_time,
        initial_window = initial_window,
        rolling_window = rolling_window,
        loss_function = loss_function,
        loss_parameter = loss_parameter,
        basis_knot_distance = basis_knot_distance,
        basis_knot_distance_power = basis_knot_distance_power,
        basis_deg = basis_deg,
        forget = forget,
        soft_threshold = soft_threshold,
        hard_threshold = hard_threshold,
        fixed_share = fixed_share,
        p_smooth_lambda = p_smooth_lambda,
        p_smooth_knot_distance = p_smooth_knot_distance,
        p_smooth_knot_distance_power = p_smooth_knot_distance_power,
        p_smooth_deg = p_smooth_deg,
        p_smooth_ndiff = p_smooth_ndiff,
        parametergrid_max_combinations = parametergrid_max_combinations,
        parametergrid = parametergrid,
        forget_past_performance = forget_past_performance,
        allow_quantile_crossing = allow_quantile_crossing,
        trace = trace
    )

    dimnames(model$specification$data$y) <- dimnames(y)

    if (intercept & !is.null(dimnames(experts))) {
        names <- c("Intercept", dimnames(experts)[[3]])
    } else if (!intercept & !is.null(dimnames(experts))) {
        names <- dimnames(experts)[[3]]
    } else if (intercept) {
        names <- c("Intercept", paste0("E", 1:dim(experts)[[3]]))
    } else if (!intercept) {
        names <- paste0("E", 1:dim(experts)[[3]])
    }

    dimnames(model$specification$data$experts)[[3]] <- names

    return(model)
}