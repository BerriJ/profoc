
# install.packages("splines2")

devtools::load_all()


Rcpp::sourceCpp(".package_dev/test.cpp")

reprex::reprex()

x <- 0:100 / 100
order <- 3
deg <- order - 1
n <- 3
mu <- 0.3
intercept <- TRUE

approx.equal()

knots <- make_knots2(n, mu = mu, deg = deg)

res <- splines2_periodic(x,
    knots = knots,
    deg = deg,
    intercept = intercept
)

ts.plot(
    as.matrix(res),
    col = 1:ncol(spline_p), lwd = 2
)

rowSums(res)


# inner_knots <- knots[(order + 1):(length(knots) - order)]

idx_inner <- (order + 1):(length(knots) - order)

idx_bounds <- c(order, length(knots) - order + 1)

bounds <- knots[c(order, length(knots) - order + 1)]

spline_p <- splines2::mSpline(x,
    knots = knots[idx_inner],
    degree = deg,
    Boundary.knots = bounds,
    periodic = TRUE,
    intercept = intercept
)


ts.plot(
    as.matrix(spline_p),
    col = 1:ncol(spline_p), lwd = 2
)

knots_seq <- c(
    bounds[1],
    inner_knots, bounds[2],
    inner_knots[1:deg] +
        bounds[2]
)

knots_seq[(1 - intercept) + order + 1:dim(spline_p)[2]]


w <- (knots_seq[(1 - intercept) + order + 1:dim(spline_p)[2]] -
    knots_seq[(1 - intercept) + 1:dim(spline_p)[2]]) /
    (order)



spline_p_new <- spline_p
for (i in 1:ncol(spline_p)) spline_p_new[, i] <- spline_p[, i] * w[i]

ts.plot(
    as.matrix(spline_p_new),
    col = 1:ncol(spline_p_new), lwd = 2
)

rowSums(spline_p_new)
