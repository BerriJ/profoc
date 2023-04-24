#' @template function_batch
#'
#' @description
#' `r lifecycle::badge("experimental")`
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
#' @param qw_crps Decides whether the sum of quantile scores (FALSE) or
#' the quantile weighted CRPS (TRUE) should be minimized.
#' Defaults to FALSE. Which corresponds to Berrisch & Ziel (2021)
#'
#' @param b_smooth A named list determining how the B-Spline matrices for
#' probabilistic smoothing are created. Default corresponds to no probabilistic
#' smoothing. See details.
#' @param p_smooth A named list determining how the hat matrices  for
#' probabilistic P-Spline smoothing are created. Default corresponds to
#' no smoothing. See details.
#'
#' @template param_forget
#'
#' @template param_soft_threshold
#' @template param_hard_threshold
#'
#' @template param_fixed_share
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
#'     p_smooth = list(lambda = 10)
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
                  b_smooth = list(
                      knots = length(tau),
                      mu = 0.5,
                      sigma = 1,
                      nonc = 0,
                      tailweight = 1,
                      deg = 1,
                      periodic = FALSE
                  ),
                  p_smooth = list(
                      knots = length(tau),
                      mu = 0.5,
                      sigma = 1,
                      nonc = 0,
                      tailweight = 1,
                      deg = 1,
                      ndiff = 1.5,
                      lambda = -Inf,
                      periodic = FALSE
                  ),
                  forget = 0,
                  soft_threshold = -Inf,
                  hard_threshold = -Inf,
                  fixed_share = 0,
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
        grid <- expand.grid(
            b_knots = val_or_def(b_smooth$knots, length(tau)),
            b_mu = val_or_def(b_smooth$mu, 0.5),
            b_sigma = val_or_def(b_smooth$sigma, 1),
            b_nonc = val_or_def(b_smooth$nonc, 0),
            b_tailweight = val_or_def(b_smooth$tailweight, 1),
            b_deg = val_or_def(b_smooth$deg, 1),
            b_periodic = val_or_def(b_smooth$periodic, FALSE),
            forget = val_or_def(forget, 0),
            soft_threshold = val_or_def(soft_threshold, -Inf),
            hard_threshold = val_or_def(hard_threshold, -Inf),
            fixed_share = val_or_def(fixed_share, 0),
            p_knots = val_or_def(p_smooth$knots, length(tau)),
            p_mu = val_or_def(p_smooth$mu, 0.5),
            p_sigma = val_or_def(p_smooth$sigma, 1),
            p_nonc = val_or_def(p_smooth$nonc, 0),
            p_tailweight = val_or_def(p_smooth$tailweight, 1),
            p_deg = val_or_def(p_smooth$deg, 1),
            p_ndiff = val_or_def(p_smooth$ndiff, 1),
            p_lambda = val_or_def(p_smooth$lambda, -Inf),
            p_periodic = val_or_def(p_smooth$periodic, FALSE)
        )

        parametergrid <- as.matrix(grid)
    } else if (ncol(parametergrid) != 20) {
        # TODO: Update Docs
        stop("Please provide a parametergrid with 20 columns.")
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

    knots <- vector("list", length = 2 * nrow(parametergrid))

    dim(knots) <- c(nrow(parametergrid), 2)

    for (i in seq_len(nrow(parametergrid))) {
        knots[i, 1][[1]] <- make_knots(
            parametergrid[i, "b_knots"],
            parametergrid[i, "b_mu"],
            parametergrid[i, "b_sigma"],
            parametergrid[i, "b_nonc"],
            parametergrid[i, "b_tailweight"],
            parametergrid[i, "b_deg"]
        )
        knots[i, 2][[1]] <- make_knots(
            parametergrid[i, "p_knots"],
            parametergrid[i, "p_mu"],
            parametergrid[i, "p_sigma"],
            parametergrid[i, "p_nonc"],
            parametergrid[i, "p_tailweight"],
            parametergrid[i, "p_deg"]
        )
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
        knots = knots,
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
