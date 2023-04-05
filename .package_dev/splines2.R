# install.packages("splines2")

devtools::load_all()

x <- 0:100 / 100
order <- 4
deg <- order - 1
n <- 3
mu <- 0.5
intercept <- TRUE

knots <- make_knots2(n, mu = mu, deg = deg)

res <- splines2_basis(x,
    knots = knots,
    deg = deg,
    periodic = TRUE,
    intercept = intercept
)

ts.plot(
    as.matrix(res),
    col = 1:ncol(res), lwd = 2
)
