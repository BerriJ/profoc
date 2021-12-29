#' @method print cppclock
#' @importFrom stats sd
#' @export
print.cppclock <- function(x, ...) {
    df <- data.frame(
        Microsecs = round(tapply(x$Nanoseconds * 1e-3, x$Name, mean)),
        SD = round(tapply(x$Nanoseconds * 1e-3, x$Name, stats::sd)),
        Obs = tapply(x$Nanoseconds, x$Name, length)
    )
    print(df, row.names = TRUE)
    invisible(df)
}