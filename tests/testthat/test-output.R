# skip_if(debug_mode)

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

model <- online(
    y = matrix(y),
    experts = experts,
    trace = FALSE
)

# Return Type
expect_type(model, "list")

# Dimensions
expect_true(all(dim(model$weights) == c(length(y) + 1, P, N)))
expect_true(all(dim(model$predictions) == c(length(y), P)))

# Missing values
expect_true(all(!is.na(model$weights)))
expect_true(all(!is.na(model$experts_loss)))
expect_true(all(!is.na(model$predictions)))
expect_true(all(!is.na(model$past_perf_wrt_params)))
expect_true(all(!is.na(model$opt_index)))
expect_true(all(!is.na(model$forecaster_loss)))
expect_true(all(!is.na(model$experts_loss)))