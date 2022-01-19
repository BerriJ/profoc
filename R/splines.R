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
            n = params[i, "n"],
            a = params[i, "beta_a"],
            b = params[i, "beta_b"],
            deg = params[i, "deg"],
            outer = params[i, "outer"]
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

make_hat_mats <- function(x,
                          n = length(x),
                          beta_a = 1,
                          beta_b = 1,
                          deg = 1,
                          outer = TRUE,
                          diff = 1.5,
                          lambda = -Inf) {
    params <- expand.grid(
        n = n,
        beta_a = beta_a,
        beta_b = beta_b,
        deg = deg,
        outer = outer,
        diff = diff,
        lambda = lambda
    )
    params <- as.matrix(params)

    hat_list <- list()

    for (i in seq_len(nrow(params))) {
        knots <- make_knots2(
            n = params[i, "n"],
            a = params[i, "beta_a"],
            b = params[i, "beta_b"],
            deg = params[i, "deg"],
            outer = params[i, "outer"]
        )
        if (params[i, "lambda"] != -Inf) {
            hat_list[[i]] <- make_hat_matrix2(
                x = x,
                knots = knots,
                deg = params[i, "deg"],
                bdiff = params[i, "diff"],
                lambda = params[i, "lambda"]
            )
        } else {
            hat_list[[i]] <- Matrix::sparseMatrix(
                i = seq_along(x),
                j = seq_along(x),
                x = 1
            )
        }
    }
    # Important for passing to C++ as arma::field object
    dim(hat_list) <- c(length(hat_list), 1)
    out <- list(hat = hat_list, params = params)
    return(out)
}