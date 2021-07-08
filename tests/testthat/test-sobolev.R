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

# We expect different results when using the sobolev space:
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = 50,
    p_smooth_ndiff = 2,
    trace = FALSE
)

boa_smooth_sobol <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = 50,
    p_smooth_ndiff = 1.5,
    trace = FALSE
)

expect_false(mean(boa_smooth_sobol$forecaster_loss) ==
    mean(boa_smooth$forecaster_loss))