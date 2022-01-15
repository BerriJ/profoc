make_basis_mats <- function(knot_distance,
                            knot_distance_power,
                            deg,
                            P_or_D) {
    sp_basis <- 1:P_or_D / (P_or_D + 1)
    params <- expand.grid(
        knot_distance = knot_distance,
        knot_distance_power = knot_distance_power,
        deg = deg
    )
    params <- as.matrix(params)
    basis_list <- list()

    for (i in seq_len(nrow(params))) {
        basis_list[[i]] <-
            make_basis_matrix(
                sp_basis,
                params[i, "knot_distance"],
                params[i, "deg"],
                params[i, "knot_distance_power"],
                P_or_D %% 2 == 0
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