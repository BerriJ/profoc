
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