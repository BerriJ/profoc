#' @importFrom stats pbeta
make_knots2 <- function(n, a, b, deg, outer = TRUE) {
    new_n <- n + ifelse(outer, 2 * (deg + 1), 0)
    if (new_n == 1) {
        seq1 <- 0.5
    } else {
        seq1 <- seq(from = 0, to = 1, length.out = new_n)
    }
    av <- mean(seq1)
    seq_beta <- pbeta(seq1, a, b)
    seq2 <- (seq_beta - av) * ((n + 2 * deg + 1) / (n + 1))
    knots <- seq2 + av
    return(knots)
}

make_basis_mats <- function(x, # Splines basis
                            n = length(x), # (vec of) Number of knots
                            beta_a = 1, # (vec of) Beta dist. alpha
                            beta_b = 1, # (vec of) Beta dist. beta
                            deg = 1, # (vec of) Degree of splines
                            outer = TRUE # (vec of) Whether to add outer knots
) {
    params <- expand.grid(
        n = n,
        beta_a = beta_a,
        beta_b = beta_b,
        deg = deg,
        outer = outer
    )
    params <- as.matrix(params)
    basis_list <- list()

    for (i in seq_len(nrow(params))) {
        knots <- make_knots2(
            n = params[i, 1],
            a = params[i, 2],
            b = params[i, 3],
            deg = params[i, 4],
            outer = params[i, 5]
        )

        basis_list[[i]] <- make_basis_matrix2(
            x = x,
            knots = knots,
            deg = params[i, 4]
        )
    }
    # Important for passing to C++ as arma::field object
    dim(basis_list) <- c(length(basis_list), 1)
    out <- list(basis = basis_list, params = params)
    return(out)
}

make_hat_mats <- function(knot_distance,
                          knot_distance_power,
                          deg,
                          lambda,
                          diff,
                          P_or_D) {
    sp_basis <- 1:P_or_D / (P_or_D + 1)
    params <- expand.grid(
        knot_distance = knot_distance,
        knot_distance_power = knot_distance_power,
        deg = deg,
        lambda = lambda,
        diff = diff
    )
    params <- as.matrix(params)

    hat_list <- list()

    for (i in seq_len(nrow(params))) {
        knots <- make_knots(
            params[i, "knot_distance"],
            params[i, "knot_distance_power"],
            params[i, "deg"],
            P_or_D %% 2 == 0
        )
        if (params[i, 4] != -Inf) {
            hat_list[[i]] <- make_hat_matrix(
                sp_basis,
                params[i, "knot_distance"],
                params[i, "lambda"],
                params[i, "diff"],
                params[i, "deg"],
                params[i, "knot_distance_power"],
                P_or_D %% 2 == 0
            )
        } else {
            hat_list[[i]] <- Matrix::sparseMatrix(
                i = 1:P_or_D,
                j = 1:P_or_D,
                x = 1
            )
        }
    }
    # Important for passing to C++ as arma::field object
    dim(hat_list) <- c(length(hat_list), 1)
    out <- list(hat = hat_list, params = params)
    return(out)
}