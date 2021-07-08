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

# We expect approximately uniform weights (concerning the quantiles), when
# lambda -> infinity
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    p_smooth_lambda = 2^40,
    trace = FALSE
)

expect_true(
    all(
        abs(
            apply(boa_smooth$weights[, , 1], 1, diff)
        ) < 1e-9
    )
)