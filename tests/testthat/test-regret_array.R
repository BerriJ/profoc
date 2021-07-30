# skip_if(debug_mode)

set.seed(1)

T <- 4
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

model_lg_true <- online(
    y = matrix(y),
    experts = experts,
    loss_gradient = TRUE,
    trace = FALSE
)

model_lg_false <- online(
    y = matrix(y),
    experts = experts,
    loss_gradient = FALSE,
    trace = FALSE
)

regret <- sweep(
    x = -model_lg_false$experts_loss,
    MARGIN = 1:2,
    FUN = "+",
    model_lg_false$forecaster_loss
)

model2 <- online(
    y = matrix(y),
    experts = experts,
    regret = regret,
    trace = FALSE
)

expect_true(
    identical(model_lg_false$weights, model2$weights)
)

expect_false(
    identical(model_lg_true$weights, model2$weights)
)

model3 <- online(
    y = matrix(y),
    experts = experts,
    regret = list(regret = regret, share = 1),
    trace = FALSE
)

expect_true(
    identical(model2$weights, model3$weights)
)

model4 <- online(
    y = matrix(y),
    experts = experts,
    regret = list(regret = regret, share = 0),
    trace = FALSE
)

expect_true(
    identical(model_lg_true$weights, model4$weights)
)

model5 <- online(
    y = matrix(y),
    experts = experts,
    regret = list(regret = regret, share = c(0, 0.5, 1)),
    trace = FALSE
)

dim(model5$past_performance)

expect_true(all(model5$past_performance[, , 1] ==
    model_lg_true$past_performance[, , 1]))

expect_true(
    all(model5$past_performance[, , 3]
    == model_lg_false$past_performance[, , 1])
)
