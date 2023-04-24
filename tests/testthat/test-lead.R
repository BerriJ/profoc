skip_if(debug_mode)
# %% Test lead for online
set.seed(1)

# Experts
N <- 2
# Observations
T <- 500
# Size of probability grid
P <- 1
prob_grid <- 1:P / (P + 1)

# Realized observations
y <- rnorm(T)

# Expert predictions
experts <- array(dim = c(T, P, N))

change_at <- 300

experts[1:(change_at - 1), , 1] <- 5
experts[1:(change_at - 1), , 2] <- -3
experts[change_at:T, , 1] <- 500
experts[change_at:T, , 2] <- -3

theor_opt2 <- (qnorm(prob_grid) - experts[T, , 2]) /
    (experts[T, , 1] - experts[T, , 2])

init_win <- 5
lead <- 10

results <- online(
    matrix(y),
    experts,
    prob_grid,
    lead_time = lead,
    trace = FALSE
)

# ts.plot(results$weights[, , 1], col = rainbow(ncol(results$weights)))
# abline(v = change_at)

expect_gt(results$weights[change_at + 5, , , 1], 0.2)

expect_lt(results$weights[change_at + lead + 5, , , 1], 0.2)
# %%

# %% Test lead for batch
results <- batch(
    matrix(y),
    experts,
    prob_grid,
    affine = TRUE,
    positive = TRUE,
    initial_window = init_win,
    lead_time = lead,
    trace = FALSE
)

# ts.plot(results$weights[, , 1], col = rainbow(ncol(results$weights)))
# abline(v = change_at)

expect_gt(results$weights[change_at + 5, , 1], 0.2)

expect_lt(results$weights[change_at + lead + 5, , 1], 0.2)
# %%
