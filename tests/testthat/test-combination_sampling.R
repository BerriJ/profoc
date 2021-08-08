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
expect_warning(boa <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    basis_knot_distance = c(0.01, 0.05, 0.5, 1),
    basis_deg = c(1, 2),
    p_smooth_lambda = c(1:2),
    parametergrid_max_combinations = 10,
    trace = FALSE
), "Warning: Too many parameter combinations possible. 10 combinations were randomly sampled. Results may depend on sampling.")

expect_true(nrow(boa$parametergrid) == 10)
