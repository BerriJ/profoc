skip_if(debug_mode)

suppressPackageStartupMessages(library(gamlss.dist))

set.seed(1)

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

# Mean
boa_mean <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    loss_function = "expectile",
    trace = FALSE
)

expect_equal(round(mean(boa_mean$predictions), 1), mean_y)

# Median
boa_median <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    trace = FALSE
)

median_y <- round(
    qSST(0.5, mu = mean_y, tau = tau_y, sigma = sd_y, nu = nu_y), 1
)

expect_equal(round(mean(boa_median$predictions), 1), median_y)