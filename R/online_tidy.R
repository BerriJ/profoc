#' Tidy the Weights of an Online object
#'
#' `tidy` will tansform the weights array of an online object
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
    return(x_tb)
}
