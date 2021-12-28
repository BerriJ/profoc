#' @method print cppclock
#' @importFrom stats sd
#' @export
print.cppclock <- function(x, ...) {
    df <- data.frame(
        Av_Nanosecs = tapply(x$Nanoseconds, x$Name, mean),
        SD = tapply(x$Nanoseconds, x$Name, stats::sd),
        Obs = tapply(x$Nanoseconds, x$Name, length)
    )
    print(df, row.names = TRUE)
    invisible(df)
}