#' Tidy the Weights of an Online object
#'
#' `tidy` will transform the weights array of an online object
#' into a tibble that is better suited for plotting and analysis.
#' @param x The weights of an `online` object.
#' @param ... Not currently used.
#' @return A tibble with columns `t` `d` `p` `k` and `w` corresponding
#' to the time, marginals, probabilities, experts, and weights
#' of the online-learning computation.
#' @importFrom tibble as_tibble
#' @importFrom generics tidy
#' @export
tidy.online.weights <- function(x, ...) {
    x_tb <- array2DF(x, responseName = "w") |> as_tibble()
    x_tb$t <- as.integer(x_tb$t)
    x_tb$p <- as.numeric(x_tb$p)
    x_tb_sorted <- x_tb[order(x_tb$t, x_tb$d, x_tb$p, x_tb$k), ]
    return(x_tb_sorted)
}

#' Tidy the Predictions of an Online object
#'
#' `tidy` will transform the `predictions` array of an online object
#' into a tibble that is better suited for plotting and analysis.
#' @param x The predictions of an `online` object.
#' @param ... Not currently used.
#' @return A tibble with columns `t` `d` `p` `k` and `w` corresponding
#' to the time, marginals, probabilities, and predictions
#' of the online-learning computation.
#' @importFrom tibble as_tibble
#' @importFrom generics tidy
#' @export
tidy.online.predictions <- function(x, ...) {
    x_tb <- array2DF(x, responseName = "prediction") |>
        as_tibble()
    x_tb$t <- as.integer(x_tb$t)
    x_tb$p <- as.numeric(x_tb$p)
    x_tb_sorted <- x_tb[order(x_tb$t, x_tb$d, x_tb$p), ]
    return(x_tb_sorted)
}

#' Tidy the Experts' losses of an Online object
#'
#' `tidy` will transform the `experts_loss` array of an online object
#' into a tibble that is better suited for plotting and analysis.
#' @param x The experts_loss of an `online` object.
#' @param ... Not currently used.
#' @return A tibble with columns `t` `d` `p` `k` and `w` corresponding
#' to the time, marginals, probabilities, and experts_loss
#' of the online-learning computation.
#' @importFrom tibble as_tibble
#' @importFrom generics tidy
#' @export
tidy.online.experts_loss <- function(x, ...) {
    x_tb <- array2DF(x, responseName = "loss") |>
        as_tibble()
    x_tb$t <- as.integer(x_tb$t)
    x_tb$p <- as.numeric(x_tb$p)
    x_tb_sorted <- x_tb[order(x_tb$t, x_tb$d, x_tb$p, x_tb$k), ]
    return(x_tb_sorted)
}

#' Tidy the Experts' losses of an Online object
#'
#' `tidy` will transform the `forecaster_loss`` array of an online object
#' into a tibble that is better suited for plotting and analysis.
#' @param x The forecaster_loss of an `online` object.
#' @param ... Not currently used.
#' @return A tibble with columns `t` `d` `p` `k` and `w` corresponding
#' to the time, marginals, probabilities, and forecaster_loss
#' of the online-learning computation.
#' @importFrom tibble as_tibble
#' @importFrom generics tidy
#' @export
tidy.online.forecaster_loss <- function(x, ...) {
    x_tb <- array2DF(x, responseName = "loss") |>
        as_tibble()
    x_tb$t <- as.integer(x_tb$t)
    x_tb$p <- as.numeric(x_tb$p)
    x_tb_sorted <- x_tb[order(x_tb$t, x_tb$d, x_tb$p), ]
    return(x_tb_sorted)
}
