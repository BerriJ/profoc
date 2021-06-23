skip_if(debug_mode)
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

model <- online(
    y = matrix(y),
    experts = experts,
    knot_distance = 1,
    trace = FALSE
)

diff(model$weights[500, , 1])

diffs <- apply(model$weights[, , ], MARGIN = 1, FUN = diff)

expect_true(sum(diffs) == 0)