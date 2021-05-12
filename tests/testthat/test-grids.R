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

# We expect that grids do effect the performance:
boa_smooth <- profoc(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    lambda = c(10, 1000),
    ndiff = c(1, 2),
    deg = c(2, 3),
    knot_distance = c(0.001, 0.01, 0.1),
    knot_distance_power = c(0.5, 1, 2),
    trace = FALSE
)

# We expect weights to sum to 1 despite the smoothing:
expect_true(all(round(apply(boa_smooth$weights, 1:2, sum), 12) == 1))

expect_true(
    all(!duplicated(apply(boa_smooth$past_perf_wrt_params, 3, mean)))
)

# We expect that the some grids do *not* affect the performance if smoothing
# is off (lambda = -Inf):
boa_pseudo_smooth <- profoc(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    ndiff = c(1, 2),
    deg = c(2, 3),
    knot_distance = c(0.001, 0.01, 0.1),
    knot_distance_power = c(0.5, 1, 2),
    trace = FALSE
)

expect_true(
    all(duplicated(apply(boa_pseudo_smooth$past_perf_wrt_params, 3, mean))[-1])
)