pgbeta <- function(x, mu = .5, sig = 1, nonc = 0) {
    if (nonc >= 0) {
        a <- mu * sig * 2
        b <- (1 - mu) * sig * 2
        c <- nonc
        return(pbeta(x, a, b, c))
    } else {
        a <- (1 - mu) * sig * 2
        b <- mu * sig * 2
        c <- abs(nonc)
        return(pbeta(1 - x, a, b, c))
    }
}

#' @importFrom stats pbeta
#' @importFrom utils head
#' @importFrom utils tail
make_knots2 <- function(n, mu = .5, sig = 1, nonc = 0, tailw = 1, deg = 1) {
    if (n < 0) {
        return(NULL)
    }
    seq_n <- seq(from = 0, to = 1, length.out = n + 2)
    seq_beta <- pgbeta(seq_n, mu, sig, nonc)
    knots <- c(
        abs(tailw) * diff(head(seq_beta, 2)) * (-deg:-1),
        seq_beta,
        1 + diff(tail(seq_beta, 2)) * (1:(deg) * abs(tailw))
    )
    return(knots)
}

make_basis_mats <- function(x, # Splines basis
                            n = length(x), # (vec of) Number of knots
                            mu = 0.5, # (vec of) Beta dist. mu
                            sigma = 1, # (vec of) Beta dist. variance param
                            nonc = 0, # (vec of) Beta dist. noncentrality
                            tailw = 1, # (vec of) Tailweight
                            deg = 1, # (vec of) Degree of splines
                            idx = NULL,
                            params = NULL) {
    if (is.null(params)) {
        params <- expand_grid_sample(
            list(
                n = n,
                mu = mu,
                sigma = sigma,
                nonc = nonc,
                tailw = tailw,
                deg = deg
            ),
            idx = idx
        )
    }

    basis_list <- list()

    for (i in seq_len(nrow(params))) {
        deg_ <- params[i, "deg"]
        n_ <- params[i, "n"]

        # Calculate temporary number of knots to enshure unity in the tails
        n_ <- n_ - (deg_ + 1) + 2
        if (n_ < 0) {
            if (deg_ + n_ >= 1) {
                print("n too small reduce deg_ by n-(deg_+1)")
                deg_ <- deg_ + n_
                n_ <- 0
            } else {
                n_ <- -1
            }
        }

        knots <- make_knots2(
            n = n_,
            mu = params[i, "mu"],
            sig = params[i, "sigma"],
            nonc = params[i, "nonc"],
            tailw = params[i, "tailw"],
            deg = deg_
        )

        if (is.null(knots)) {
            basis_list[[i]] <- Matrix::Matrix(rep.int(1, length(x)), sparse = TRUE)
        } else {
            basis_list[[i]] <- make_basis_matrix2(
                x = x,
                knots = knots,
                deg = deg_
            )
        }
    }
    # Important for passing to C++ as arma::field object
    dim(basis_list) <- c(length(basis_list), 1)
    out <- list(basis = basis_list, params = params)
    return(out)
}

make_hat_mats <- function(x,
                          n = length(x),
                          mu = 0.5,
                          sigma = 1,
                          nonc = 0,
                          tailw = 1,
                          deg = 1,
                          ndiff = 1.5,
                          lambda = -Inf,
                          idx = NULL,
                          params = NULL) {
    if (is.null(params)) {
        params <- expand_grid_sample(
            list(
                n = n,
                mu = mu,
                sigma = sigma,
                nonc = nonc,
                tailw = tailw,
                deg = deg,
                ndiff = ndiff,
                lambda = lambda
            ),
            idx = idx
        )
    }

    hat_list <- list()

    for (i in seq_len(nrow(params))) {
        deg_ <- params[i, "deg"]
        n_ <- params[i, "n"]

        # Calculate temporary number of knots to enshure unity in the tails
        n_ <- n_ - (deg_ + 1) + 2

        if (n_ < 0) {
            if (deg_ + n_ >= 1) {
                print("n too small reduce deg_ by n-(deg_+1)")
                deg_ <- deg_ + n_
                n_ <- 0
            } else {
                print("n=1, constant case")
                n_ <- -1
            }
        }

        knots <- make_knots2(
            n = n_,
            mu = params[i, "mu"],
            sig = params[i, "sigma"],
            nonc = params[i, "nonc"],
            tailw = params[i, "tailw"],
            deg = deg_
        )

        if (params[i, "lambda"] != -Inf) {
            hat_list[[i]] <- make_hat_matrix2(
                x = x,
                knots = knots,
                deg = deg_,
                bdiff = params[i, "ndiff"],
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
