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

# %% Online setting

grids <- list()

grids$general <- as.matrix(expand.grid(
    "forget_regret" = 0,
    "soft_threshold" = 0,
    "hard_threshold" = 0,
    "fixed_share" = 0,
    "basis_pr_idx" = c(1, 2),
    "basis_mv_idx" = 1,
    "hat_pr_idx" = 1,
    "hat_mv_idx" = 1,
    "gamma" = 1,
    "loss_share" = 0,
    "regret_share" = 0
))

grids$b_smooth_pr <- as.matrix(expand.grid(
    n = 9,
    mu = c(0.1, 0.2),
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1
))

grids$b_smooth_mv <- as.matrix(expand.grid(
    n = 9,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1
))

grids$p_smooth_pr <- as.matrix(expand.grid(
    n = 1,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1,
    diff = 1.5,
    lambda = -Inf
))

grids$p_smooth_mv <- as.matrix(expand.grid(
    n = 1,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1,
    diff = 1.5,
    lambda = -Inf
))

mod <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    parametergrids = grids,
    trace = FALSE
)

expect_true(all(mod$parametergrid[, "basis_pr_idx"] == grids$general[, "basis_pr_idx"]))

# %%
