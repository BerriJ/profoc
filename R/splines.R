make_basis_mats <- function(knot_distance,
                            knot_distance_power,
                            deg,
                            P_or_D) {
    sp_basis <- 1:P_or_D / (P_or_D + 1)
    params <- expand.grid(knot_distance, knot_distance_power, deg)
    basis_list <- list()

    for (i in seq_len(nrow(params))) {
        knots <- make_knots(
            params[i, 1],
            params[i, 2],
            params[i, 3],
            P_or_D %% 2 == 0
        )
        basis_list[[i]] <-
            make_basis_matrix2(sp_basis, knots, params[i, 3])
    }
    # Important for passing to C++ as arma::field object
    dim(basis_list) <- c(length(basis_list), 1)
    return(basis_list)
}