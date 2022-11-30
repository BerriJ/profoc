
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
    model$weights <- list_to_array(model$weights)
    model$past_performance <- list_to_array(model$past_performance)
    model$experts_loss <- list_to_array(model$experts_loss)
    dimnames(model$experts_loss)[[4]] <- names$experts[[4]]
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
        warning(
            paste(
                "Warning: Too many parameter combinations possible.",
                n,
                "combinations were randomly sampled. Results depend on sampling."
            )
        )
    }

    if (is.null(idx)) {
        idx <- sort(sample.int(N, n, replace = FALSE))
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
