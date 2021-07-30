skip_if(debug_mode)

# %% Test Setup
set.seed(1)

# Experts
N <- 2
# Observations
T <- 50
# Size of probability grid
P <- 9
prob_grid <- 1:P / (P + 1)

threshold_val <- 0.3

dev <- c(-1, 4)
experts_sd <- c(1, 1)

# Realized observations
y <- rnorm(n = T)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = dev[1], sd = experts_sd[1])
    experts[t, , 2] <- qnorm(prob_grid, mean = dev[2], sd = experts_sd[2])
}
# %%

# %% Create grid

grid <- expand.grid(
    0.01, # basis_kstep
    1, # kstep_power
    1, # basis_deg
    0.1, # forget_regret
    c(0.01, 0.15), # soft_threshold
    0, # hard_threshold
    c(0, 0.04), # fixed_share
    -Inf, # lambda
    0.01, # p_smooth_knot_distance
    1, # p_smooth_knot_distance_power
    3, # p_smooth_deg
    2, # p_smooth_ndiff
    1, # gamma
    0, # loss_share
    0 # regret_share
)
grid <- as.matrix(grid)
# %%

# %% Check if custom grid is used
res <- online(
    matrix(y),
    experts,
    parametergrid = grid,
    trace = FALSE
)

expect_true(
    all(res$parametergrid == grid)
)
# %%

# %% Check dimension error
expect_error(online(
    matrix(y),
    experts,
    parametergrid = grid[, -1],
    trace = FALSE
), "Please provide a parametergrid with 15 columns.")
# %%

# %% Batch setting
res <- batch(
    matrix(y),
    experts,
    parametergrid = grid[, 1:12], # No gamma parameter in batch
    trace = FALSE
)

expect_true(
    all(res$parametergrid == grid[, 1:12])
)
# %%

# %% Check dimension error
expect_error(batch(
    matrix(y),
    experts,
    parametergrid = grid[, -c(12:13)], # No gamma parameter in batch
    trace = FALSE
), "Please provide a parametergrid with 12 columns.")
# %%