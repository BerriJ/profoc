
# install.packages("splines2")

x <- 0:100 / 100
order <- 2
deg <- order - 1
knots <- c(0.00, 0.082, 0.23, 0.47, 1.00)

spline_p <- splines2::mSpline(x,
    knots = knots[c(-1, -5)],
    degree = deg,
    Boundary.knots = c(0, 1),
    periodic = TRUE,
    intercept = TRUE
)

ts.plot(
    as.matrix(spline_p),
    col = 1:ncol(spline_p), lwd = 2
)

knots_simple <- c(0, knots, 1 + knots[2])

w <- c()
for (i in 2:5) {
    w[i - 1] <- (knots_simple[deg + 1 + i] - knots_simple[i]) / (deg + 1)
}

spline_p_new <- spline_p
for (i in 1:ncol(spline_p)) spline_p_new[, i] <- spline_p[, i] * w[i]

ts.plot(
    as.matrix(spline_p_new),
    col = 1:ncol(spline_p_new), lwd = 2
)

rowSums(spline_p_new)
