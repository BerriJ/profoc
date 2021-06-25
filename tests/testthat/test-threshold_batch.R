skip_if(debug_mode)

# %% Test threshold batch - unconstrained
set.seed(1)

# Experts
N <- 2
# Observations
T <- 100
# Size of probability grid
P <- 1
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

results_unconstrained <- batch(
    matrix(y),
    experts,
    prob_grid,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_unconstrained_hard <- batch(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_unconstrained_soft <- batch(
    matrix(y),
    experts,
    prob_grid,
    soft_threshold = threshold_val,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

# ts.plot(results_unconstrained$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Results_unconstrained as is"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_unconstrained_hard$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Hard"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_unconstrained_soft$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Test threshold batch - convex
results_convex <- batch(
    matrix(y),
    experts,
    prob_grid,
    positive = TRUE,
    affine = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_convex_hard <- batch(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    positive = TRUE,
    affine = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_convex_soft <- batch(
    matrix(y),
    experts,
    prob_grid,
    soft_threshold = threshold_val,
    positive = TRUE,
    affine = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

# ts.plot(results_convex$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Results_convex as is"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_convex_hard$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Hard Convex"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_convex_soft$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft Convex"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Test threshold batch - convex + intercept
results_convex_intercept <- batch(
    matrix(y),
    experts,
    prob_grid,
    positive = TRUE,
    affine = TRUE,
    intercept = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_convex_intercept_hard <- batch(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    positive = TRUE,
    affine = TRUE,
    intercept = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

results_convex_intercept_soft <- batch(
    matrix(y),
    experts,
    prob_grid,
    soft_threshold = threshold_val,
    positive = TRUE,
    affine = TRUE,
    intercept = TRUE,
    trace = FALSE,
    basis_knot_distance = 0.01,
    basis_deg = 1
)

# ts.plot(results_convex_intercept$weights[, 1, ],
#     col = 1:3,
#     ylim = c(0, 1),
#     main = "Results convex + intercept as is"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_convex_intercept_hard$weights[, 1, ],
#     col = 1:3,
#     ylim = c(0, 1),
#     main = "Threshold Hard Convex + Intercept"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_convex_intercept_soft$weights[, 1, ],
#     ylim = c(0, 1),
#     col = 1:3,
#     main = "Threshold Soft Convex + Intercept"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Test intercept influnence in convex setting

# Test if intercept weights are in the threshold region
expect_true(
    any(results_convex_intercept$weights[, , 1] < threshold_val &
        results_convex_intercept$weights[, , 1] > 0.01)
)

# Thresholds shouldn't influence intercept
expect_true(
    any(threshold_val < results_convex_intercept_hard$weights[, , 1] &
        results_convex_intercept_hard$weights[, , 1] > 0.01)
)

# %%

# %% Thresholds should influence experts in convex setting
# Test if expert weights are in the threshold region
expect_true(
    any(results_convex_intercept$weights[, , 3] < threshold_val &
        results_convex_intercept$weights[, , 3] > 0.01)
)

# Weights removed from the threshold region? (either set to > threshold or 0)
expect_false(
    any(results_convex_intercept_hard$weights[, , 3] < threshold_val &
        results_convex_intercept_hard$weights[, , 3] > 0.01)
)
# %%