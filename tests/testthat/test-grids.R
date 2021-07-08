skip_if(debug_mode)

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

# We expect that grids do effect the performance:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = c(10, 1000),
    p_smooth_ndiff = c(1, 2),
    p_smooth_deg = c(2, 3),
    p_smooth_knot_distance = c(0.001, 0.1),
    p_smooth_knot_distance_power = c(0.5, 1, 2),
    basis_deg = 1,
    basis_knot_distance = 0.01,
    basis_knot_distance_power = 1,
    trace = FALSE
)

# We expect weights to sum to 1 despite the smoothing:
expect_true(all(round(apply(boa_smooth$weights, 1:2, sum), 12) == 1))

expect_true(
    all(!duplicated(apply(boa_smooth$past_performance, 3, mean)))
)

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    basis_knot_distance = 0.1,
    p_smooth_lambda = c(1),
    p_smooth_ndiff = seq(from = 1, to = 2, by = 0.2),
    trace = FALSE
)

expect_true(
    all(!duplicated(apply(boa_smooth$past_performance, 3, mean)))
)