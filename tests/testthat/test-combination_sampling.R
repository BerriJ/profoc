skip_if(debug_mode)

set.seed(1)

T <- 32
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

expect_message(boa <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        knots = c(15, 20, 30, 50, 99),
        deg = c(1, 2),
        lambda = c(1:2)
    ),
    parametergrid_max_combinations = 3,
    trace = FALSE
), "Too many parameter combinations possible. 3 combinations were randomly sampled. Results depend on sampling.")

# Test if the seed works as expected
set.seed(1)
suppressMessages({
    first <- online(
        y = matrix(y),
        tau = prob_grid,
        experts = experts,
        soft_threshold = c(0, 0.1, 0.2, 0.3, 0.4),
        forget_regret = c(0, 0.1, 0.2, 0.3, 0.4),
        parametergrid_max_combinations = 3,
        trace = FALSE
    )
})

set.seed(1)
suppressMessages({
    second <- online(
        y = matrix(y),
        tau = prob_grid,
        experts = experts,
        soft_threshold = c(0, 0.1, 0.2, 0.3, 0.4),
        forget_regret = c(0, 0.1, 0.2, 0.3, 0.4),
        parametergrid_max_combinations = 3,
        trace = FALSE
    )
})

expect_true(all(first$parametergrid == second$parametergrid))

set.seed(1)
suppressMessages({
    first <- online(
        y = matrix(y),
        tau = prob_grid,
        experts = experts,
        soft_threshold = c(0, 0.1, 0.2, 0.3, 0.4),
        forget_regret = c(0, 0.1, 0.2, 0.3, 0.4),
        parametergrid_max_combinations = 3,
        trace = FALSE
    )
})

set.seed(2)
suppressMessages({
    second <- online(
        y = matrix(y),
        tau = prob_grid,
        experts = experts,
        soft_threshold = c(0, 0.1, 0.2, 0.3, 0.4),
        forget_regret = c(0, 0.1, 0.2, 0.3, 0.4),
        parametergrid_max_combinations = 3,
        trace = FALSE
    )
})

expect_false(all(first$parametergrid == second$parametergrid))
