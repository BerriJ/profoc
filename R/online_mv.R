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
#' @param smooth_pr tbd
#' @param smooth_mv tbd
#'
#' @param forget_regret Share of past regret not to be considered, resp. to be
#' forgotten in every iteration of the algorithm. Defaults to 0.
#'
#' @template param_soft_threshold
#' @template param_hard_threshold
#'
#' @template param_fixed_share
#'
#'
#' @param gamma Scaling parameter for the learning rate.
#'
#' @template param_parametergrid_max_combinations
#' @template param_parametergrid
#' @template param_forget_past_performance
#'
#' @template param_allow_quantile_crossing
#'
#' @param init A named list containing "init_weights": Matrix of dimension 1xK or PxK used as starting weights. 1xK represents the constant solution with equal weights over all P, whereas specifying a PxK matrix allows different starting weights for each P. "R0" a matrix of dimension PxK or 1xK used as starting regret.
#' @param loss User specified loss array. Can also be a list with elements "loss_array"
#' and "share", share mixes the provided loss with the loss calculated by
#' profoc. 1 means, only the provided loss will be used. share can also be
#' vector of shares to consider.
#' @param regret User specified regret array. If specific, the regret will
#'  not be calculated by profoc. Can also be a list with elements "regret_array"
#' and "share", share mixes the provided regret with the regret calculated by
#' profoc. 1 means, only the provided regret will be used. share can also be
#' vector of shares to consider.
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
#' }
#' @importFrom abind asub
#' @export
online_mv <- function(y, experts, tau,
                      lead_time = 0,
                      loss_function = "quantile",
                      loss_parameter = 1,
                      loss_gradient = TRUE,
                      method = "bewa",
                      smooth_pr = list(
                          basis = list(
                              knot_distance = 1 / (P + 1),
                              knot_distance_power = 1,
                              deg = 1
                          ),
                          penalized = list(
                              lambda = -Inf,
                              knot_distance = 1 / (P + 1),
                              knot_distance_power = 1,
                              deg = 1,
                              ndiff = 1.5
                          )
                      ),
                      smooth_mv = list(
                          basis = list(
                              knot_distance = 1 / (D + 1),
                              knot_distance_power = 1,
                              deg = 1
                          ),
                          penalized = list(
                              lambda = -Inf,
                              knot_distance = 1 / (D + 1),
                              knot_distance_power = 1,
                              deg = 1,
                              ndiff = 1.5
                          )
                      ),
                      forget_regret = 0,
                      soft_threshold = -Inf,
                      hard_threshold = -Inf,
                      fixed_share = 0,
                      gamma = 1,
                      parametergrid_max_combinations = 100,
                      parametergrid = NULL,
                      forget_past_performance = 0,
                      allow_quantile_crossing = FALSE,
                      init = NULL,
                      loss = NULL,
                      regret = NULL,
                      trace = TRUE) {
    edim <- dim(experts)

    if (is.vector(y)) {
        y <- matrix(y)
    }

    if (length(edim) == 3) {
        if (ncol(y) > 1) { # multivariate point
            experts <- array(unlist(experts), dim = c(edim[1], edim[2], 1, edim[3]))
            experts <- lapply(seq_len(edim[1]),
                asub,
                x = experts,
                dims = 1,
                drop = TRUE
            )
            dim(experts) <- c(edim[1], 1)
        } else if (ncol(y) == 1) { # Probabilistic univariate
            experts <- lapply(seq_len(edim[1]),
                asub,
                x = experts,
                dims = 1,
                drop = FALSE
            )
            dim(experts) <- c(edim[1], 1)
        }
    } else if (length(edim) == 4) { # multivariate probabilistic
        experts <- lapply(seq_len(edim[1]),
            asub,
            x = experts,
            dims = 1,
            drop = TRUE
        )
        dim(experts) <- c(edim[1], 1)
    }
    exdim <- dim(experts[[1]])

    T <- dim(experts)[1]
    D <- dim(experts[[1]])[1]
    P <- dim(experts[[1]])[2]
    K <- dim(experts[[1]])[3]

    if (nrow(experts) - nrow(y) < 0) {
        stop("Number of provided expert predictions has to match or exceed observations.")
    }

    if (nrow(y) <= lead_time) {
        stop("Number of expert predictions need to exceed lead_time.")
    }

    if (is.null(loss)) {
        loss_array <- array(, dim = c(0, 0, 0))
        loss_share <- 0
    } else if (is.array(loss)) {
        loss_array <- loss
        loss_share <- 1
    } else if (is.list(loss)) {
        loss_array <- loss$loss
        loss_share <- loss$share
    }

    if (is.null(regret)) {
        regret_array <- array(, dim = c(0, 0, 0))
        regret_share <- 0
    } else if (is.array(regret)) {
        regret_array <- regret
        regret_share <- 1
    } else if (is.list(regret)) {
        regret_array <- regret$regret
        regret_share <- regret$share
    }

    val_or_def <- function(val, def) {
        if (is.null(val)) {
            return(def)
        } else {
            return(val)
        }
    }

    if (is.null(parametergrid)) {
        grid <- expand.grid(
            val_or_def(
                smooth_pr$basis$knot_distance,
                1 / (P + 1)
            ),
            val_or_def(
                smooth_pr$basis$knot_distance_power,
                1
            ),
            val_or_def(
                smooth_pr$basis$deg,
                1
            ),
            forget_regret,
            soft_threshold,
            hard_threshold,
            fixed_share,
            val_or_def(
                smooth_pr$penalized$lambda,
                -Inf
            ),
            val_or_def(
                smooth_pr$penalized$knot_distance,
                1 / (P + 1)
            ),
            val_or_def(
                smooth_pr$penalized$knot_distance_power,
                1
            ),
            val_or_def(
                smooth_pr$penalized$deg,
                1
            ),
            val_or_def(
                smooth_pr$penalized$ndiff,
                1.5
            ),
            gamma,
            loss_share,
            regret_share,
            val_or_def(
                smooth_mv$basis$knot_distance,
                1 / (P + 1)
            ),
            val_or_def(
                smooth_mv$basis$knot_distance_power,
                1
            ),
            val_or_def(
                smooth_mv$basis$deg,
                1
            ),
            val_or_def(
                smooth_mv$penalized$lambda,
                -Inf
            ),
            val_or_def(
                smooth_mv$penalized$knot_distance,
                1 / (P + 1)
            ),
            val_or_def(
                smooth_mv$penalized$knot_distance_power,
                1
            ),
            val_or_def(
                smooth_mv$penalized$deg,
                1
            ),
            val_or_def(
                smooth_mv$penalized$ndiff,
                1.5
            )
        )

        parametergrid <- as.matrix(grid)
    } else if (ncol(parametergrid) != 15) {
        stop("Please provide a parametergrid with 15 columns.")
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

    if (is.null(init$init_weights)) {
        init$init_weights <- matrix(
            1 / exdim[3],
            nrow = exdim[2],
            ncol = exdim[3]
        )
    } else if (nrow(init$init_weights) == 1) {
        init$init_weights <- matrix(init$init_weights,
            nrow = exdim[2],
            ncol = exdim[3],
            byrow = TRUE
        )
    } else if (
        (nrow(init$init_weights) != 1 &
            nrow(init$init_weights) != exdim[2]) |
            ncol(init$init_weights) != exdim[3]) {
        stop("Either a 1xK or PxK matrix of initial weights must be supplied.")
    }
    init$init_weights <- pmax(init$init_weights, exp(-350))
    init$init_weights <- init$init_weights / rowSums(init$init_weights)

    if (is.null(init$R0)) {
        init$R0 <- matrix(0,
            nrow = exdim[2],
            ncol = exdim[3],
        )
    } else if (nrow(init$R0) == 1) {
        init$R0 <- matrix(init$R0,
            nrow = exdim[2],
            ncol = exdim[3],
            byrow = TRUE
        )
    } else if (
        (nrow(init$R0) != 1 &
            nrow(init$R0) != exdim[2]) |
            ncol(init$R0) != exdim[3]) {
        stop("R0 must be 1xK or PxK.")
    }

    model <- online_rcpp_mv(
        y = y,
        experts = experts,
        tau = tau,
        lead_time = lead_time,
        loss_function = loss_function,
        loss_parameter = loss_parameter,
        loss_gradient = loss_gradient,
        method = method,
        param_grid = parametergrid,
        forget_past_performance = forget_past_performance,
        allow_quantile_crossing = allow_quantile_crossing,
        w0 = init$init_weights,
        R0 = init$R0,
        loss_array = loss_array,
        regret_array = regret_array,
        trace = trace
    )

    dimnames(model$specification$data$y) <- dimnames(y)

    tmp <- array(
        NA,
        c(
            dim(model$weights)[1],
            dim(model$weights[[1]])[1],
            dim(model$weights[[1]])[2],
            dim(model$weights[[1]])[3]
        )
    )
    for (i in seq_len(dim(model$weights)[1])) tmp[i, , , ] <- model$weights[[i]]
    model$weights <- tmp

    tmp <- array(
        NA,
        c(
            dim(model$past)[1],
            dim(model$past[[1]])[1],
            dim(model$past[[1]])[2],
            dim(model$past[[1]])[3]
        )
    )
    for (i in seq_len(dim(model$past)[1])) tmp[i, , , ] <- model$past[[i]]
    model$past <- tmp

    tmp <- array(
        NA,
        c(
            dim(model$experts_loss)[1],
            dim(model$experts_loss[[1]])[1],
            dim(model$experts_loss[[1]])[2],
            dim(model$experts_loss[[1]])[3]
        )
    )
    for (i in seq_len(dim(model$experts_loss)[1])) tmp[i, , , ] <- model$experts_loss[[i]]
    model$experts_loss <- tmp

    parnames <- c(
        "basis_knot_distance",
        "basis_knot_distance_power",
        "basis_deg",
        "forget_regret",
        "threshold_soft",
        "threshold_hard",
        "fixed_share",
        "p_smooth_lambda",
        "p_smooth_knot_distance",
        "p_smooth_knot_distance_power",
        "p_smooth_deg",
        "smooth_diff",
        "gamma",
        "loss_share",
        "regret_share",
        "mv_basis_knot_distance",
        "mv_basis_knot_distance_power",
        "mv_basis_deg",
        "mv_p_smooth_lambda",
        "mv_p_smooth_knot_distance",
        "mv_p_smooth_knot_distance_power",
        "mv_p_smooth_deg",
        "mv_p_smooth_ndiff"
    )

    colnames(model$chosen_parameters) <- parnames
    colnames(model$parametergrid) <- parnames

    return(model)
}