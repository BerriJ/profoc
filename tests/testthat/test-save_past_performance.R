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

# We expect that past_performance is reported:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        lambda = c(10, 1000),
        ndiff = c(1, 2),
        deg = c(2, 3)
    ),
    save_past_performance = TRUE,
    trace = FALSE
)

expect_true(all(dim(boa_smooth$past_performance) == c(1000, 1, 99, 8)))

expect_true(
    all(round(boa_smooth$past_performance[T, 1, 50, ], 7) ==
        c(
            0.2728130, 0.3878642, 0.3796495, 0.3793674,
            0.3309897, 0.3310782, 0.2906660, 0.2907555
        ))
)

expect_true(boa_smooth$specification$parameters$save_past_performance)

# The default of save_past_performance is expected to be FALSE
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

expect_true(all(dim(boa_smooth$past_performance) == c(1, 0, 0, 0)))

expect_false(boa_smooth$specification$parameters$save_past_performance)
