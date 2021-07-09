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

boa_smooth_crossed <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    sorting = "none",
    trace = FALSE
)

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    sorting = "ex_post",
    trace = FALSE
)

# We expect that quantile crossing happens and therefore losses are different:
expect_true(
    sum(boa_smooth_crossed$forecaster_loss) != sum(boa_smooth$forecaster_loss)
)

# We expect that the "sorted loss" is less than the "unsorted loss"
expect_lt(
    sum(boa_smooth$forecaster_loss), sum(boa_smooth_crossed$forecaster_loss)
)

# Multivariate setting (sorting will be disabled automatically)
T <- 1000
N <- 2
P <- 2
prob_grid <- c(0.5, 0.5)

mean_y <- 0
sd_y <- 5

# Realized observations
y <- matrix(c(rnorm(n = T), rnorm(n = T)), ncol = 2)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = -5, sd = 2)
    experts[t, , 2] <- qnorm(prob_grid, mean = 5, sd = 2)
}

boa_smooth_crossed <- online(
    y = y,
    tau = prob_grid,
    experts = experts,
    sorting = "none",
    trace = FALSE
)

expect_warning(
    boa_smooth <- online(
        y = y,
        tau = prob_grid,
        experts = experts,
        sorting = "ex_post",
        trace = FALSE
    )
)


# We expect that the "sorted loss" is less than the "unsorted loss"
expect_equal(
    sum(boa_smooth$forecaster_loss), sum(boa_smooth_crossed$forecaster_loss)
)