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
        lambda = c(100, 1000),
        ndiff = c(1, 2),
        deg = c(2, 3)
    ),
    save_predictions_grid = TRUE,
    trace = FALSE
)

expect_true(all(dim(boa_smooth$specification$objects$predictions_grid) == c(1000, 1, 1)))

osize_full <- object.size(boa_smooth$specification$objects$predictions_grid)

smpl <- boa_smooth$specification$objects$predictions_grid[[563]][, 50, ]

expect_equal(
    smpl,
    c(
        -0.16726231, -0.16427737, -0.21685071, -0.22053253,
        -0.19778184, -0.19708188, -0.19235012, -0.20075103
    ),
    tolerance = 0.0000001
)

expect_true(boa_smooth$specification$parameters$save_predictions_grid)

# The default of save_predictions_grid is expected to be FALSE

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_pr = list(
        lambda = c(100, 1000),
        ndiff = c(1, 2),
        deg = c(2, 3)
    ),
    trace = FALSE
)

expect_true(all(dim(boa_smooth$specification$objects$predictions_grid) == c(1, 1, 1)))

osize_small <- object.size(boa_smooth$specification$objects$predictions_grid)

smpl <- boa_smooth$specification$objects$predictions_grid[[1]][, 50, ]

expect_equal(
    smpl,
    c(
        -0.07082431, -0.07394584, -0.14018784, -0.13642743,
        -0.03533878, -0.03516170, -0.11598618, -0.11580723
    ),
    tolerance = 0.000001
)


expect_false(boa_smooth$specification$parameters$save_predictions_grid)

# The predictions grid has to expand if lead_time >= 1 even if save_predictions_grid = FALSE.

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    lead_time = 3,
    experts = experts,
    p_smooth_pr = list(
        lambda = c(10, 1000),
        ndiff = c(1, 2),
        deg = c(2, 3)
    ),
    trace = FALSE
)

expect_true(all(dim(boa_smooth$specification$objects$predictions_grid) == c(4, 1, 1)))

osize_lead <- object.size(boa_smooth$specification$objects$predictions_grid)

expect_false(boa_smooth$specification$parameters$save_predictions_grid)

expect_true(osize_full > osize_lead)
expect_true(osize_lead > osize_small)
