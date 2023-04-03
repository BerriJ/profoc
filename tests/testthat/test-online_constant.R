skip_if(debug_mode)

set.seed(1)

T <- 1000
N <- 2
P <- 99
prob_grid <- 1:P / (P + 1)

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
    tau = prob_grid,
    b_smooth_pr = list(knots = -1),
    trace = FALSE
)

diff(model$weights[500, , , 1])

diffs <- apply(model$weights[, , , ], MARGIN = 1, FUN = diff)

expect_true(sum(diffs) == 0)

expect_true(dim(model$specification$objects$basis_pr[[1]])[2] == 1)
