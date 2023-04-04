# Simple setting of hat function

devtools::load_all()
devtools::load_all()

x <- 1:19 / 20
lambda <- 10
order <- 2
deg <- order - 1
n_inner <- 10 # Inner knots
mu <- 0.5

knots <- make_knots2(n_inner, deg = deg)
image(as.matrix(make_hat_matrix2(x, knots, deg, 1, lambda = 1, periodic = TRUE)))


par(mfrow = c(3, 5))
for (i in 2:4) {
    for (m in c(0.2, 0.4, 0.5, 0.6, 0.8)) {
        order <- i
        deg <- order - 1
        mu <- m
        knots <- make_knots2(n_inner, mu = mu, deg = deg)
        # ts.plot(
        #     splines2_periodic(x, knots, deg),
        #     col = seq_along(knots),
        #     lwd = 2
        # )
        P <- as.matrix(make_hat_matrix2(x, knots, deg, 1, lambda = 5, periodic = TRUE))
        image(P, main = paste("Deg=", deg, ", mu=", mu))
    }
}


# par(mfrow = c(1, 2))

# B <- splines2_basis(x, knots, deg)

# dim(B)

# B_p <- splines2_periodic(x, knots, deg)

# dim(B_p)

# ts.plot(
#     B_p,
#     col = seq_along(knots),
#     lwd = 2
# )

# P <- penalty_periodic(knots, order)

# image(P)

# hat <- solve(t(B_p) %*% B_p + lambda * P) %*% t(B_p)

# image(hat)

# par(mfrow = c(1, 2))

# x <- cumsum(rnorm(99))
# # x <- c(x, rev(x))

# plot(x, type = "l")

# fitted <- B_p %*% solve(t(B_p) %*% B_p + lambda * P) %*% t(B_p) %*% matrix(x)

# plot(fitted, type = "l")
