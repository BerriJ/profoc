skip_if(debug_mode)
suppressPackageStartupMessages(library(gamlss.dist))

set.seed(1)

# Experts
N <- 2
# Observations
T <- 1000
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

experts[, , 1] <- 5
experts[, , 2] <- -5


# Mean linear
o_weights_linear <- oracle(
    matrix(y),
    experts,
    loss_function = "expectile",
    tau = prob_grid,
)

expect_equal(mean(o_weights_linear$predictions), mean_y, tolerance = 0.2)

# Mean Convex
convex <- oracle(
    matrix(y),
    experts,
    loss_function = "expectile",
    tau = prob_grid,
    positive = TRUE,
    affine = TRUE
)

expect_equal(mean(convex$predictions), mean_y, tolerance = 0.2)

# Median - linear
o_weights_linear <- oracle(
    matrix(y),
    experts,
    loss_function = "quantile",
    tau = prob_grid
)

median_y <- qSST(0.5, mu = mean_y, tau = tau_y, sigma = sd_y, nu = nu_y)

expect_equal(mean(o_weights_linear$predictions), median_y, tolerance = 0.2)

# Median - convex
o_weights_convex <- oracle(
    matrix(y),
    experts,
    loss_function = "quantile",
    tau = prob_grid,
    positive = TRUE,
    affine = TRUE
)

expect_equal(mean(o_weights_convex$predictions), median_y, tolerance = 0.2)
