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

model1 <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    method = "boa",
    trace = FALSE
)

# Return Type
expect_type(model1, "list")

# Dimensions
expect_true(all(dim(model1$weights) == c(length(y) + 1, 1, P, N)))
expect_true(all(dim(model1$predictions) == c(length(y), 1, P)))

# Missing values
expect_true(all(!is.na(model1$weights)))
expect_true(all(!is.na(model1$experts_loss)))
expect_true(all(!is.na(model1$predictions)))
expect_true(all(!is.na(model1$past_perf_wrt_params)))
expect_true(all(!is.na(model1$opt_index)))
expect_true(all(!is.na(model1$forecaster_loss)))
expect_true(all(!is.na(model1$experts_loss)))

expect_true(round(mean(model1$forecaster_loss), 10) == 0.6763033692)

expect_true(
    all(
        round(tail(model1$forecaster_loss)[, , 80], 10) ==
            c(0.7660281621, 0.8528983718, 0.9663739408, 0.8740971730, 1.0867901304, 0.9297034346)
    )
)

model2 <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    method = "bewa",
    trace = FALSE
)

# ombination 0.6806536

# Return Type
expect_type(model2, "list")

# Dimensions
expect_true(all(dim(model2$weights) == c(length(y) + 1, 1, P, N)))
expect_true(all(dim(model2$predictions) == c(length(y), 1, P)))

# Missing values
expect_true(all(!is.na(model2$weights)))
expect_true(all(!is.na(model2$experts_loss)))
expect_true(all(!is.na(model2$predictions)))
expect_true(all(!is.na(model2$past_perf_wrt_params)))
expect_true(all(!is.na(model2$opt_index)))
expect_true(all(!is.na(model2$forecaster_loss)))
expect_true(all(!is.na(model2$experts_loss)))

expect_true(round(mean(model2$forecaster_loss), 10) == 0.6806535578)

expect_true(
    all(
        round(tail(model2$forecaster_loss)[, , 80], 10) ==
            c(0.7660281622, 0.8528983719, 0.9663739409, 0.8740971731, 1.0867901305, 0.9297034347)
    )
)

model3 <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    method = "ewa",
    gamma = 0.1,
    trace = FALSE
)

model3$specification$parameters$method

# Return Type
expect_type(model3, "list")

# Dimensions
expect_true(all(dim(model3$weights) == c(length(y) + 1, 1, P, N)))
expect_true(all(dim(model3$predictions) == c(length(y), 1, P)))

# Missing values
expect_true(all(!is.na(model3$weights)))
expect_true(all(!is.na(model3$experts_loss)))
expect_true(all(!is.na(model3$predictions)))
expect_true(all(!is.na(model3$past_perf_wrt_params)))
expect_true(all(!is.na(model3$opt_index)))
expect_true(all(!is.na(model3$forecaster_loss)))
expect_true(all(!is.na(model3$experts_loss)))

expect_true(round(mean(model3$forecaster_loss), 10) == 0.7067959411)

expect_true(
    all(
        round(tail(model3$forecaster_loss)[, , 80], 10) ==
            c(0.7660281623, 0.8528983720, 0.9663739409, 0.8740971732, 1.0867901305, 0.9297034348)
    )
)

model4 <- online(
    y = matrix(y),
    experts = experts,
    tau = prob_grid,
    method = "ml_poly",
    trace = FALSE
)

# Return Type
expect_type(model4, "list")

# Dimensions
expect_true(all(dim(model4$weights) == c(length(y) + 1, 1, P, N)))
expect_true(all(dim(model4$predictions) == c(length(y), 1, P)))

# Missing values
expect_true(all(!is.na(model4$weights)))
expect_true(all(!is.na(model4$experts_loss)))
expect_true(all(!is.na(model4$predictions)))
expect_true(all(!is.na(model4$past_perf_wrt_params)))
expect_true(all(!is.na(model4$opt_index)))
expect_true(all(!is.na(model4$forecaster_loss)))
expect_true(all(!is.na(model4$experts_loss)))

expect_true(round(mean(model4$forecaster_loss), 10) == 0.7135789929)

expect_true(
    all(
        round(tail(model4$forecaster_loss)[, , 80], 10) ==
            c(0.7660281623, 0.8528983720, 0.9663739409, 0.8740971732, 1.0867901305, 0.9297034348)
    )
)
