#' @template function_online
#'
#' @template param_y
#' @template param_experts
#' @template param_tau
#'
#' @template param_lead_time
#'
#' @template param_loss_function
#' @template param_loss_parameter
#' @param loss_gradient Determines if a linearized version of the loss is used.
#'
#' @template param_method
#'
#' @template param_basis_knot_distance_online
#' @template param_basis_knot_distance_power
#' @template param_basis_deg_online
#'
#' @param forget_regret Share of past regret not to be considered, resp. to be
#' forgotten in every iteration of the algorithm. Defaults to 0.
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
#' @param gamma Scaling parameter for the learning rate.
#'
#' @template param_parametergrid_max_combinations
#' @template param_parametergrid
#' @template param_forget_past_performance
#'
#' @template param_allow_quantile_crossing
#'
#' @param init_weights Matrix of dimension 1xK or PxK used as starting weights. 1xK represents the constant solution with equal weights over all P, whereas specifying a PxK matrix allows different starting weights for each P.
#' @param loss_array User specified loss array. If specified, the loss will not be calculated by profoc.
#' @param regret_array User specified regret array. If specific, the regret will not be calculated by profoc.
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
#' model <- online(
#'     y = matrix(y),
#'     experts = experts,
#'     p_smooth_lambda = 10
#' )
#'
#' print(model)
#' plot(model)
#' autoplot(model)
#'
#' new_y <- matrix(rnorm(1)) # Realized
#' new_experts <- experts[T, , , drop = FALSE]
#'
#' # Update will update the model object, no need for new assignment
#' update(model, new_y = new_y, new_experts = new_experts)
#'
#' # Use predict to combine new_experts, model$predictions will be extended
#' predict(model, new_experts = new_experts)
#' @export
online <- function(y, experts,
                   tau = 1:dim(experts)[2] / (dim(experts)[2] + 1),
                   lead_time = 0,
                   loss_function = "quantile",
                   loss_parameter = 1,
                   loss_gradient = TRUE,
                   method = "bewa",
                   basis_knot_distance = 1 / (dim(experts)[2] + 1),
                   basis_knot_distance_power = 1,
                   basis_deg = 1,
                   forget_regret = 0,
                   soft_threshold = -Inf,
                   hard_threshold = -Inf,
                   fixed_share = 0,
                   p_smooth_lambda = -Inf,
                   p_smooth_knot_distance = basis_knot_distance,
                   p_smooth_knot_distance_power = basis_knot_distance_power, p_smooth_deg = basis_deg,
                   p_smooth_ndiff = 1.5,
                   gamma = 1,
                   parametergrid_max_combinations = 100,
                   parametergrid = NULL,
                   forget_past_performance = 0,
                   allow_quantile_crossing = FALSE,
                   init_weights = NULL,
                   loss_array = NULL,
                   regret_array = NULL,
                   trace = TRUE) {



    # Ensure that online_rcpp does not expand a grid for basis_knot_distance
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

    if (is.null(loss_array)) {
        loss_array <- array(, dim = c(0, 0, 0))
    }

    if (is.null(regret_array)) {
        regret_array <- array(, dim = c(0, 0, 0))
    }

    model <- online_rcpp(
        y = y, experts = experts, tau = tau,
        lead_time = lead_time,
        loss_function = loss_function,
        loss_parameter = loss_parameter,
        loss_gradient = loss_gradient,
        method = method,
        basis_knot_distance = basis_knot_distance,
        basis_knot_distance_power = basis_knot_distance_power,
        basis_deg = basis_deg,
        forget_regret = forget_regret,
        soft_threshold = soft_threshold,
        hard_threshold = hard_threshold,
        fixed_share = fixed_share,
        p_smooth_lambda = p_smooth_lambda,
        p_smooth_knot_distance = p_smooth_knot_distance,
        p_smooth_knot_distance_power = p_smooth_knot_distance_power,
        p_smooth_deg = p_smooth_deg,
        p_smooth_ndiff = p_smooth_ndiff,
        gamma = gamma,
        parametergrid_max_combinations = parametergrid_max_combinations,
        parametergrid = parametergrid,
        forget_past_performance = forget_past_performance,
        allow_quantile_crossing = allow_quantile_crossing,
        init_weights = init_weights,
        loss_array = loss_array,
        regret_array = regret_array,
        trace = trace
    )

    dimnames(model$specification$data$y) <- dimnames(y)

    dimnames(model$specification$data$experts) <- dimnames(experts)

    if (is.null(dimnames(model$specification$data$experts)[[3]])) {
        dimnames(model$specification$data$experts)[[3]] <-
            paste0("E", 1:dim(experts)[[3]])
    }

    return(model)
}