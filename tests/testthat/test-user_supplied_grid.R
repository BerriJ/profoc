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
    b_knots = length(prob_grid),
    b_mu = 0.5,
    b_sigma = 1,
    b_nonc = 0,
    b_tailweight = 1,
    b_deg = 1,
    b_periodic = FALSE,
    forget = 0,
    soft_threshold = -Inf,
    hard_threshold = -Inf,
    fixed_share = 0,
    p_knots = length(prob_grid),
    p_mu = 0.5,
    p_sigma = 1,
    p_nonc = 0,
    p_tailweight = 1,
    p_deg = c(1, 2, 3),
    p_ndiff = c(1, 2),
    p_lambda = c(-Inf, 1, 5, 10),
    p_periodic = c(TRUE, FALSE)
)

grid <- as.matrix(grid)
# %%

# %% Batch setting
res <- batch(
    matrix(y),
    experts,
    parametergrid = grid, # No gamma parameter in batch
    trace = FALSE
)
# %%

# %% Check dimension error
expect_error(batch(
    matrix(y),
    experts,
    parametergrid = grid[, -12], # No gamma parameter in batch
    trace = FALSE
), "Please provide a parametergrid with 20 columns.")
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
    deg = 1,
    periodic = FALSE
))

grids$b_smooth_mv <- as.matrix(expand.grid(
    n = 9,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1,
    periodic = FALSE
))

grids$p_smooth_pr <- as.matrix(expand.grid(
    n = 1,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1,
    diff = 1.5,
    lambda = -Inf,
    periodic = FALSE
))

grids$p_smooth_mv <- as.matrix(expand.grid(
    n = 1,
    mu = 0.5,
    sigma = 1,
    nonc = 0,
    tailw = 1,
    deg = 1,
    diff = 1.5,
    lambda = -Inf,
    periodic = FALSE
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
