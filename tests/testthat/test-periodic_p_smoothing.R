set.seed(1)

T <- 1000
N <- 2
P <- 99
prob_grid <- 1:P / (P + 1)

mean_y <- 0
sd_y <- 5

# Realized observations
y <- rnorm(n = T)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = -5, sd = 2)
    experts[t, , 2] <- qnorm(prob_grid, mean = 5, sd = 2)
}

# We expect that grids do affects the performance:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        knots = 5,
        lambda = 1,
        ndiff = 1,
        deg = 2,
        periodic = FALSE
    ),
    trace = FALSE
)

dn <- boa_smooth$weights[T + 1, 1, 1, 2] - boa_smooth$weights[T + 1, 1, P, 2]

expect_equal(round(dn, 6), 0.362049)

# We expect that grids do affects the performance:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        knots = 5,
        lambda = 1,
        ndiff = 1,
        deg = 2,
        periodic = TRUE
    ),
    trace = FALSE
)

# plot(prob_grid, boa_smooth$weights[T + 1, 1, , 2], type = "o", lwd = 2)

dp <- boa_smooth$weights[T + 1, 1, 1, 2] - boa_smooth$weights[T + 1, 1, P, 2]

expect_equal(round(dp, 6), 0.028608)
