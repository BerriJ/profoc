suppressPackageStartupMessages(library(gamlss.dist))

# Experts
N <- 2
# Observations
T <- 10000
# Size of probability grid
P <- 1
prob_grid <- 1:P / (P + 1)

mean_y <- 0
sd_y <- 5
tau_y <- 6
nu_y <- 2

# Realized observations
y <- rSST(n = T, mu = mean_y, tau = tau_y, sigma = sd_y, nu = nu_y)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = -5, sd = 2)
    experts[t, , 2] <- qnorm(prob_grid, mean = 5, sd = 2)
}

# Mean linear
o_weights_linear <- oracle(
    matrix(y),
    experts,
    loss_function = "expectile",
    tau = prob_grid,
    convex_constraint = FALSE
)

expect_equal(mean(o_weights_linear$predictions), mean_y, tolerance = 0.1)

# Mean Convex
convex <- oracle(
    matrix(y),
    experts,
    loss_function = "expectile",
    tau = prob_grid,
    convex_constraint = TRUE
)

expect_equal(mean(convex$predictions), mean_y, tolerance = 0.1)

# Median - linear
o_weights_linear <- oracle(
    matrix(y),
    experts,
    loss_function = "quantile",
    tau = prob_grid,
    convex_constraint = FALSE
)

median_y <- qSST(0.5, mu = mean_y, tau = tau_y, sigma = sd_y, nu = nu_y)

expect_equal(mean(o_weights_linear$predictions), median_y, tolerance = 0.1)

# Median - convex
o_weights_convex <- oracle(
    matrix(y),
    experts,
    loss_function = "quantile",
    tau = prob_grid,
    convex_constraint = TRUE
)

expect_equal(mean(o_weights_convex$predictions), median_y, tolerance = 0.1)