#' Print method for online models
#'
#' Prints the average loss of all experts and the forecast combination.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @export
print.online <- function(x, ...) {
    print_common(x)
}