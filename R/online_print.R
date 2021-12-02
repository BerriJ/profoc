#' Print method for online models
#'
#' Prints the average loss of all experts and the forecast combination.
#' @param x Object of class inheriting from 'online'
#' @param ...  further arguments are ignored
#' @export
print.online <- function(x, ...) {
    expert_names <- dimnames(x$experts_loss)[[4]]
    experts_loss <- round(apply(x$experts_loss, 4, mean), 5)
    forecaster_loss <- mean(x$forecaster_loss)

    print(data.frame(
        Name = c(expert_names, "Combination"),
        Loss = c(experts_loss, forecaster_loss)
    ), row.names = FALSE, right = FALSE)
}