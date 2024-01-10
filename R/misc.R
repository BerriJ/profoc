get_seed <- function() {
    sample.int(.Machine$integer.max, 1)
}

#' Create experts list to be used in conline class
#'
#' This function works in conjunction with the conline class.
#' It takes a matrix of experts and a matrix of outcomes and
#' returns a list of experts which fulfills all properties
#' that are needed for passing it to the an instance of conline.
#' @param experts array of predictions with dimension T x D x P x K
#' (Observations x Variables x Quantiles x Experts) or T x D x K or T x P x K.
#' @param y  A matrix of outcomes with dimension T x D.
#' @param output_with_names Defaults to FALSE. If TRUE, the function
#' returns a list with
#' the experts list, the names of the variables (dnames) and the
#' names of the experts (enames).
#' @importFrom abind asub adrop
#' @export
init_experts_list <- function(experts, y, output_with_names = FALSE) {
    edim <- dim(experts)
    if (length(edim) == 3) {
        enames <- dimnames(experts)[[3]]
        if (is.null(enames)) {
            enames <- paste0("E", 1:edim[3])
        }
        if (ncol(y) > 1) { # multivariate point
            if (is.null(dimnames(experts)[[2]])) {
                dnames <- paste0("D", 1:edim[2])
            } else {
                dnames <- dimnames(experts)[[2]]
            }
            experts <- array(
                unlist(experts),
                dim = c(edim[1], edim[2], 1, edim[3])
            )
            experts <- lapply(seq_len(edim[1]),
                asub,
                x = experts,
                dims = 1,
                drop = FALSE
            )
            experts <- lapply(experts, adrop, drop = 1)
            dim(experts) <- c(edim[1], 1)
        } else if (ncol(y) == 1) { # univariate probabilistic
            dnames <- "D1"
            experts <- lapply(seq_len(edim[1]),
                asub,
                x = experts,
                dims = 1,
                drop = FALSE
            )
            dim(experts) <- c(edim[1], 1)
        }
    } else if (length(edim) == 4) { # multivariate probabilistic
        if (is.null(dimnames(experts)[[2]])) {
            dnames <- paste0("D", 1:edim[2])
        } else {
            dnames <- dimnames(experts)[[2]]
        }
        enames <- dimnames(experts)[[4]]
        if (is.null(enames)) {
            enames <- paste0("E", 1:edim[4])
        }
        experts <- array_to_list(experts)
    }

    if (output_with_names) {
        return(list(experts = experts, dnames = dnames, enames = enames))
    } else {
        return(experts)
    }
}

#' @importFrom abind asub adrop
array_to_list <- function(x) {
    x <- lapply(seq_len(dim(x)[1]),
        function(y, x) adrop(asub(x, idx = y, dims = 1, FALSE), drop = 1),
        x = x
    )
    dim(x) <- c(length(x), 1, 1)
    return(x)
}

list_to_array <- function(x) {
    dx <- dim(x)
    dx1 <- dim(x[[1]])
    out <- array(NA, c(dx[1], dx1[1], dx1[2], dx1[3]))
    for (i in seq_len(dx[1])) out[i, , , ] <- x[[i]]
    return(out)
}

darken <- function(color, factor = 1.4) {
    col <- col2rgb(color)
    col <- col / factor
    col <- rgb(t(col), maxColorValue = 255)
    col
}

val_or_def <- function(val, def) {
    if (is.null(val)) {
        return(def)
    } else {
        return(val)
    }
}


#' Post Process Data from conline Class
#'
#' This function works in conjunction with the conline class.
#' After the main learning task, it takes the output of the
#' conline class and returns an object suitable for, visualization,
#' further, and deployment.
#' analysis.
#' @param model_instance  An instance of conline.
#' @param names A named list with dimnames of `y` and `experts`.
#' @export
post_process_model <- function(model_instance, names) {
    # Generate output
    model <- list(
        predictions = model_instance$predictions,
        predictions_got_sorted = model_instance$predictions_got_sorted,
        weights = model_instance$weights,
        forecaster_loss = model_instance$loss_for,
        experts_loss = model_instance$loss_exp,
        past_performance = model_instance$past_performance,
        opt_index = model_instance$opt_index + 1, # Respect one-based indexing
        parametergrid = model_instance$params,
        params_basis_pr = model_instance$params_basis_pr,
        params_basis_mv = model_instance$params_basis_mv,
        params_hat_pr = model_instance$params_hat_pr,
        params_hat_mv = model_instance$params_hat_mv
    )

    model[["specification"]] <-
        list(
            data =
                list(
                    y = model_instance$y,
                    experts = model_instance$experts,
                    tau = model_instance$tau
                ),
            objects =
                list(
                    weights_tmp = model_instance$weights_tmp,
                    predictions_grid = model_instance$predictions_grid,
                    cum_performance = model_instance$cum_performance,
                    hat_pr = model_instance$hat_pr,
                    hat_mv = model_instance$hat_mv,
                    basis_pr = model_instance$basis_pr,
                    basis_mv = model_instance$basis_mv,
                    V = model_instance$V,
                    E = model_instance$E,
                    eta = model_instance$eta,
                    R = model_instance$R,
                    beta = model_instance$beta,
                    beta0field = model_instance$beta0field
                ),
            parameters =
                list(
                    lead_time = model_instance$lead_time,
                    loss_function = model_instance$loss_function,
                    loss_parameter = model_instance$loss_parameter,
                    loss_gradient = model_instance$loss_gradient,
                    method = model_instance$method,
                    forget_past_performance = model_instance$forget_past_performance,
                    allow_quantile_crossing = model_instance$allow_quantile_crossing,
                    save_past_performance = model_instance$save_past_performance,
                    save_predictions_grid = model_instance$save_predictions_grid
                )
        )

    attr(model, "class") <- c("online", "list")

    model$specification$data$names <- names
    dimnames(model$specification$data$y) <- names$y
    model$past_performance <- list_to_array(model$past_performance)
    model$experts_loss <- list_to_array(model$experts_loss)
    dimnames(model$experts_loss)[[3]] <- model$specification$data$tau
    dimnames(model$experts_loss)[[4]] <- names$experts[[4]]

    # Post process weights
    model$weights <- list_to_array(model$weights)

    dimnames(model$weights) <- list(
        t = 1:dim(model$weights)[1],
        d = names$experts[[2]],
        p = names$experts[[3]],
        k = names$experts[[4]]
    )

    class(model$weights) <- "online.weights"

    dimnames(model$predictions) <- list(
        t = 1:dim(model$predictions)[1],
        d = names$experts[[2]],
        p = names$experts[[3]]
    )

    class(model$predictions) <- "online.predictions"

    dimnames(model$forecaster_loss) <- list(
        t = 1:dim(model$forecaster_loss)[1],
        d = names$experts[[2]],
        p = names$experts[[3]]
    )

    class(model$forecaster_loss) <- "online.forecaster_loss"

    dimnames(model$experts_loss) <- list(
        t = 1:dim(model$experts_loss)[1],
        d = names$experts[[2]],
        p = names$experts[[3]],
        k = names$experts[[4]]
    )

    class(model$experts_loss) <- "online.experts_loss"

    return(model)
}

expand_grid_sample <- function(vecs,
                               n = NULL,
                               idx = NULL,
                               verbose = FALSE) {
    v_len <- sapply(vecs, length)
    N <- prod(v_len)

    if (is.null(idx) & is.null(n)) {
        return(as.matrix(expand.grid(vecs)))
    } else if (!is.null(idx) & !is.null(n)) {
        stop("Only one of n or idx can be specified")
    } else if (is.null(n)) {
        n <- length(idx)
    }

    if (n >= N) {
        return(as.matrix(expand.grid(vecs)))
    } else if (verbose == TRUE) {
        message(
            paste(
                "Too many parameter combinations possible.",
                n,
                "combinations were randomly sampled. Results depend on sampling."
            )
        )
    }

    if (is.null(idx)) {
        if (N <= 10^15) {
            idx <- sort(sample.int(N, n, replace = FALSE))
        } else {
            idx <- sample_int(N, n, get_seed())
        }
    }

    grid_sample <- sapply(
        seq_along(vecs),
        function(i) {
            vecs[[i]][((ceiling(idx / prod(v_len[seq_len(i - 1)])) - 1)
            %% (v_len[i]) + 1)]
        }
    )
    colnames(grid_sample) <- names(vecs)
    return(grid_sample)
}
