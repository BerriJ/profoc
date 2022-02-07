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
#' @param qw_crps Decides wether the sum of quantile scores (FALSE) or
#' the quantile weighted CRPS (TRUE) should be minimized.
#' Defaults to FALSE. Which corresponds to Berrisch & Ziel (2021)
#'
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
#' \dontrun{
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
#' }
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
                  qw_crps = FALSE,
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
    if (nrow(y) <= initial_window) {
        stop("Initial estimation window greater or equal to input data.")
    }

    if (initial_window > rolling_window) {
        stop("Initial estimation window bigger than rolling_window.")
    }

    if (nrow(experts) - nrow(y) < 0) {
        stop("Number of provided expert predictions has to match or exceed observations.")
    }

    if (ncol(y) > 1 & !allow_quantile_crossing) {
        warning("Warning: allow_quantile_crossing set to true since multivariate prediction target was provided.")
        # Bool is set inside C++
    }
    if (is.null(parametergrid)) {
        if (missing(p_smooth_knot_distance)) {
            p_smooth_knot_distance <- 0
            inh_kstep <- TRUE
        } else {
            inh_kstep <- FALSE
        }

        if (missing(p_smooth_knot_distance_power)) {
            p_smooth_knot_distance_power <- 0
            inh_kstep_p <- TRUE
        } else {
            inh_kstep_p <- FALSE
        }

        if (missing(p_smooth_deg)) {
            p_smooth_deg <- 0
            inh_deg <- TRUE
        } else {
            inh_deg <- FALSE
        }

        grid <- expand.grid(
            basis_knot_distance,
            basis_knot_distance_power,
            basis_deg,
            forget,
            soft_threshold,
            hard_threshold,
            fixed_share,
            p_smooth_lambda,
            p_smooth_knot_distance,
            p_smooth_knot_distance_power,
            p_smooth_deg,
            p_smooth_ndiff
        )

        if (inh_kstep) {
            grid[, 9] <- grid[, 1]
        }

        if (inh_kstep_p) {
            grid[, 10] <- grid[, 2]
        }

        if (inh_deg) {
            grid[, 11] <- grid[, 3]
        }

        parametergrid <- as.matrix(grid)
    } else if (ncol(parametergrid) != 12) {
        stop("Please provide a parametergrid with 12 columns.")
    }


    if (nrow(parametergrid) > parametergrid_max_combinations) {
        warning(
            paste(
                "Warning: Too many parameter combinations possible.",
                parametergrid_max_combinations,
                "combinations were randomly sampled. Results may depend on sampling."
            )
        )
        parametergrid <- parametergrid[sample(
            x = 1:nrow(parametergrid),
            size = parametergrid_max_combinations
        ), ]
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
        qw_crps = qw_crps,
        param_grid = parametergrid,
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