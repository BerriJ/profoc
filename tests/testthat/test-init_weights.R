# skip_if(debug_mode)

set.seed(1)

T <- 50
N <- 2
P <- 9
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

# Initil weights should be uniform on default:
model <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    trace = FALSE
)

expect_true(all(model$weights[1, , , ] == 0.5))

# Weights should be populated for all quantiles if not specified individually
init_weights <- matrix(c(0.3, 0.7), byrow = T, ncol = N, nrow = P)
init_weights <- array(init_weights, dim = c(1, P, N))

model <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    init = list(init_weights = init_weights),
    trace = FALSE
)

expect_true(all(model$weights[1, , , ] == init_weights[1, , ]))

# Weights can be specified for each quantile individually:
init_weights <- matrix(nrow = P, ncol = N)
init_weights[, 1] <- 1:9 / 10
init_weights[, 2] <- 9:1 / 10
init_weights <- array(init_weights, dim = c(1, P, N))

model <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    init = list(init_weights = init_weights),
    trace = FALSE
)

expect_true(all(model$weights[1, , , ] == init_weights[1, , ]))

# Weights should allways sum to 1:
init_weights <- matrix(nrow = P, ncol = N)
init_weights[, 1] <- 1:9 / 10
init_weights[, 2] <- 9:1 / 10

init_weights[1:5, ] <- init_weights[1:5, ] * 2
init_weights[6:9, ] <- init_weights[6:9, ] / 3
init_weights <- array(init_weights, dim = c(1, P, N))

model <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    init = list(init_weights = init_weights),
    trace = FALSE
)

expect_true(all(rowSums(model$weights[1, , , ]) == 1))