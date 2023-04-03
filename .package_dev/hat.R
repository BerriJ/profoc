# Simple setting of hat function

devtools::load_all()
devtools::load_all()

x <- 1:99 / 100
lambda <- 10
order <- 2
deg <- order - 1
n_inner <- 50 # Inner knots
mu <- 0.5

knots <- make_knots2(n_inner, deg = deg)

par(mfrow = c(1, 2))

B <- splines2_basis(x, knots, deg)

dim(B)

B_p <- splines2_periodic(x, knots, deg)

dim(B_p)

ts.plot(
    B_p,
    col = seq_along(knots),
    lwd = 2
)

P <- penalty_periodic(knots, order)

image(P)

hat <- solve(t(B_p) %*% B_p + lambda * P) %*% t(B_p)

image(hat)

par(mfrow = c(1, 2))

x <- cumsum(rnorm(99))
# x <- c(x, rev(x))

plot(x, type = "l")

fitted <- B_p %*% solve(t(B_p) %*% B_p + lambda * P) %*% t(B_p) %*% matrix(x)

plot(fitted, type = "l")

make_hat_matrix2(x, knots, deg, 1.5, lambda = 1, periodic = FALSE)
