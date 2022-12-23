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

# We expect that grids do affects the performance:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        lambda = c(10, 1000),
        ndiff = c(1, 2),
        deg = c(2, 3)
    ),
    trace = FALSE
)

# We expect weights to sum to 1 despite the smoothing:
expect_true(all(round(apply(boa_smooth$weights, 1:3, sum), 13) == 1))

expect_true(
    all(!duplicated(apply(boa_smooth$past_performance, 3, mean)))
)

# Enshure that development does not affect the performance:
expect_true(
    all(
        round(tail(boa_smooth$forecaster_loss)[, , 80], 7) ==
            c(0.1551044, 0.2390565, 0.3468112, 0.2466774, 0.4535191, 0.2875674)
    )
)

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        lambda = c(1),
        ndiff = seq(from = 1, to = 2, by = 0.2)
    ),
    trace = FALSE
)

expect_true(
    all(!duplicated(apply(boa_smooth$past_performance, 3, mean)))
)

# Test forget_past_performance
# We expect that forget affects the performance:
without_forget <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        knots = c(5, 10, 20),
        lambda = c(1, 10, 100, 1000),
        ndiff = c(1, 1.5, 2)
    ),
    forget_past_performance = 0,
    trace = FALSE
)

loss_without <- mean(without_forget$forecaster_loss)
smooth_mat_without <- as.numeric(without_forget$opt_index)

with_forget <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        knots = c(5, 10, 20),
        lambda = c(1, 10, 100, 1000),
        ndiff = c(1, 1.5, 2)
    ),
    forget_past_performance = 0.01,
    trace = FALSE
)

loss_with <- mean(with_forget$forecaster_loss)
smooth_mat_with <- as.numeric(with_forget$opt_index)

expect_true(loss_without != loss_with)
expect_true(!all(smooth_mat_without == smooth_mat_with))
