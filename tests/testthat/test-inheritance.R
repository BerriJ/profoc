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

# We expect inheritance from basis > smooth:
mod <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = c(10),
    basis_deg = c(1.3, 1.7),
    basis_knot_distance = c(0.043, 0.447),
    basis_knot_distance_power = c(0.7, 1.3),
    trace = FALSE
)

expect_true(
    identical(
        mod$parametergrid[, "basis_knot_distance"],
        mod$parametergrid[, "p_smooth_knot_distance"]
    )
)
expect_true(
    identical(
        mod$parametergrid[, "basis_deg"],
        mod$parametergrid[, "p_smooth_deg"]
    )
)
expect_true(
    identical(
        mod$parametergrid[, "basis_knot_distance_power"],
        mod$parametergrid[, "p_smooth_knot_distance_power"]
    )
)

# We do not expect inheritance from basis > smooth:
mod <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = c(10),
    p_smooth_deg = c(1.3, 1.7),
    p_smooth_knot_distance = c(0.043, 0.447),
    p_smooth_knot_distance_power = c(0.7, 1.3),
    trace = FALSE
)

expect_false(
    identical(
        mod$parametergrid[, "basis_knot_distance"],
        mod$parametergrid[, "p_smooth_knot_distance"]
    )
)
expect_false(
    identical(
        mod$parametergrid[, "basis_deg"],
        mod$parametergrid[, "p_smooth_deg"]
    )
)
expect_false(
    identical(
        mod$parametergrid[, "basis_knot_distance_power"],
        mod$parametergrid[, "p_smooth_knot_distance_power"]
    )
)