pgbeta <- function(x, mu = .5, sig = 1, nonc = 0) {
    a <- (1 - mu) * sig * 2
    b <- mu * sig * 2
    c <- abs(nonc)
    seq_beta <- pbeta(x, a, b, c)
    if (nonc < 0) seq_beta <- rev(1 - seq_beta)
    return(seq_beta)
}

#' Create a vector of knots for splines
#'
#' This function creates a knot vector for splines. The knots are distributed
#' according to a beta distribution. The first input defines the number of inner
#' knots. The total number of knots is \code{n + 2 * order}.
#'
#' @param n Number of knots
#' @param mu Beta distribution location parameter
#' @param sig Beta distribution scale parameter
#' @param nonc Beta distribution noncentrality parameter
#' @param tailw Tailweight
#' @param deg Degree of splines
#' @importFrom stats pbeta
#' @importFrom utils head
#' @importFrom utils tail
#' @export
make_knots <- function(n, mu = .5, sig = 1, nonc = 0, tailw = 1, deg = 1) {
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

#' Create a List of Basis Matrices
#'
#' This function creates a list of basis matrices and the corresponding
#' parameters. It is used in `online()` to create the basis matrices
#' for basis smoothing.
#' @param x The predictor variable
#' @param n Number of knots
#' @param mu Beta distribution location parameter
#' @param sigma Beta distribution scale parameter
#' @param nonc Beta distribution noncentrality parameter
#' @param tailw Tailweight
#' @param deg Degree of splines
#' @param periodic Create periodic basis
#' @param idx `make_basis_mats()` will create a grid containing all
#' combinations of the parameters. If idx is set, this grid will
#' be subsetted to the rows specified by idx.
#' @param params Instead of the arguments above, a grid (data.frame
#' or named matrix) of parameters can be passed directly.
#' @export
make_basis_mats <- function(x, # Splines basis
                            n = length(x), # (vec of) Number of knots
                            mu = 0.5, # (vec of) Beta dist. mu
                            sigma = 1, # (vec of) Beta dist. variance param
                            nonc = 0, # (vec of) Beta dist. noncentrality
                            tailw = 1, # (vec of) Tailweight
                            deg = 1, # (vec of) Degree of splines
                            periodic = FALSE, # Create periodic splines
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
                periodic = periodic
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
                print("n too small reduce deg by n-(deg+1)")
                deg_ <- deg_ + n_
                n_ <- 0
            } else {
                n_ <- -1
            }
        }

        knots <- make_knots(
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
            basis_list[[i]] <- make_basis_matrix(
                x = x,
                knots = knots,
                deg = deg_,
                periodic = params[i, "periodic"]
            )
        }
    }
    # Important for passing to C++ as arma::field object
    dim(basis_list) <- c(length(basis_list), 1)
    out <- list(basis = basis_list, params = params)
    return(out)
}

#' Create a List of Hat Matrices
#'
#' This function creates a list of hat matrices and the corresponding
#' parameters. It is used in `online()` to create the hat matrices
#' for penalized smoothing.
#' @param x The predictor variable
#' @param n Number of knots
#' @param mu Beta distribution location parameter
#' @param sigma Beta distribution scale parameter
#' @param nonc Beta distribution noncentrality parameter
#' @param tailw Tailweight
#' @param deg Degree of splines
#' @param ndiff Sets the degree of the differencing matrix for creating
#' the penalty
#' @param lambda Penalty parameter (higher values lead to higher penalty)
#' @param periodic Create periodic penalty
#' @param idx `make_hat_mats()` will create a grid containing all
#' combinations of the parameters. If idx is set, this grid will
#' be subsetted to the rows specified by idx.
#' @param params Instead of the arguments above, a grid (data.frame
#' or named matrix) of parameters can be passed directly.
#' @export
make_hat_mats <- function(x,
                          n = length(x),
                          mu = 0.5,
                          sigma = 1,
                          nonc = 0,
                          tailw = 1,
                          deg = 1,
                          ndiff = 1.5,
                          lambda = -Inf,
                          periodic = FALSE,
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
                lambda = lambda,
                periodic = periodic
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

        knots <- make_knots(
            n = n_,
            mu = params[i, "mu"],
            sig = params[i, "sigma"],
            nonc = params[i, "nonc"],
            tailw = params[i, "tailw"],
            deg = deg_
        )

        if (params[i, "lambda"] != -Inf) {
            hat_list[[i]] <- make_hat_matrix(
                x = x,
                knots = knots,
                deg = deg_,
                bdiff = params[i, "ndiff"],
                lambda = params[i, "lambda"],
                periodic = params[i, "periodic"]
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
