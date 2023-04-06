# Simple setting of hat function
devtools::load_all()
devtools::load_all()

# %%
x <- 1:100 / 100
lambda <- 10
order <- 4
deg <- order - 1
n_inner <- 5 # Inner knots
mu <- 0.3
knots <- make_knots2(n_inner, deg = deg, mu = mu)

# B <- splines2_basis(x, knots, deg, intercept = TRUE, periodic = TRUE)

H <- as.matrix(make_hat_matrix2(x, knots, deg, 1, lambda = 5, periodic = TRUE))

image(H)

par(mfrow = c(1, 1))
set.seed(3)
real <- cumsum(rnorm(length(x)))
plot(real, type = "l")

fitted <- H %*% matrix(real)
plot(c(fitted, fitted), type = "o")
# %%

# %% Vis hat matrices

x <- 1:19 / 20
lambda <- 5
n_inner <- 10 # Inner knots

par(mfrow = c(3, 5))
for (i in 2:4) {
    for (m in c(0.2, 0.4, 0.5, 0.6, 0.8)) {
        order <- i
        deg <- order - 1
        mu <- m
        knots <- make_knots2(n_inner, mu = mu, deg = deg)
        P <- as.matrix(make_hat_matrix2(x, knots, deg, 1, lambda = lambda, periodic = TRUE))
        image(P, main = paste("Deg=", deg, ", mu=", mu))
    }
}
# %%
