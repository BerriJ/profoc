get_seed <- function() {
    sample.int(.Machine$integer.max, 1)
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

post_process_model <- function(model, names) {
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
