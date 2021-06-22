skip_if(debug_mode)

# %% Test lead for online
set.seed(1)

# Experts
N <- 2
# Observations
T <- 50
# Size of probability grid
P <- 1
prob_grid <- 1:P / (P + 1)

threshold_val <- 0.3

dev <- c(-1, 4)
experts_sd <- c(1, 1)

# Realized observations
y <- rnorm(n = T)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = dev[1], sd = experts_sd[1])
    experts[t, , 2] <- qnorm(prob_grid, mean = dev[2], sd = experts_sd[2])
}

results <- online(
    matrix(y),
    experts,
    prob_grid,
    trace = FALSE
)

results_hard <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    trace = FALSE
)

results_soft <- online(
    matrix(y),
    experts,
    prob_grid,
    soft_threshold = threshold_val,
    trace = FALSE
)

# ts.plot(results$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Results as is"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_hard$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Hard"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_soft$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft"
# )
# abline(h = threshold_val, col = "grey")
# %%


# %% Mixture
set.seed(1)

results_soft_hard <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    soft_threshold = threshold_val,
    trace = FALSE
)

results_soft_hard_all_post <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    soft_threshold = threshold_val,
    ex_post_threshold_hard = TRUE,
    ex_post_threshold_soft = TRUE,
    trace = FALSE
)

# ts.plot(results_soft_hard$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft-Hard"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_soft_hard_all_post$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft-Hard Ex-Post"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Mix Ex-Post
set.seed(1)

results_soft_hard_hard_post <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    soft_threshold = threshold_val,
    ex_post_threshold_hard = TRUE,
    trace = FALSE
)

results_soft_hard_soft_post <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    soft_threshold = threshold_val,
    ex_post_threshold_soft = TRUE,
    trace = FALSE
)

# ts.plot(results_soft_hard_hard_post$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft-Hard with Hard Ex-Post"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_soft_hard_soft_post$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft-Hard with Soft Ex-Post"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Ex-Post setting
set.seed(1)

results_hard_post <- online(
    matrix(y),
    experts,
    prob_grid,
    hard_threshold = threshold_val,
    ex_post_threshold_hard = TRUE,
    trace = FALSE
)

results_soft_post <- online(
    matrix(y),
    experts,
    prob_grid,
    soft_threshold = threshold_val,
    ex_post_threshold_soft = TRUE,
    trace = FALSE
)

# ts.plot(results_hard_post$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Hard Ex-Post"
# )
# abline(h = threshold_val, col = "grey")

# ts.plot(results_soft_post$weights[, 1, ],
#     ylim = c(0, 1),
#     main = "Threshold Soft Ex-Post"
# )
# abline(h = threshold_val, col = "grey")
# %%

# %% Tests
losses <- c()

# Threshold results differ from normal results
losses[1] <- mean(results$forecaster_loss)
losses[2] <- mean(results_hard$forecaster_loss)
losses[3] <- mean(results_soft$forecaster_loss)
losses[4] <- mean(results_hard_post$forecaster_loss)
losses[5] <- mean(results_soft_post$forecaster_loss)
losses[6] <- mean(results_soft_hard$forecaster_loss)
losses[7] <- mean(results_soft_hard_all_post$forecaster_loss)
losses[8] <- mean(results_soft_hard_hard_post$forecaster_loss)
losses[9] <- mean(results_soft_hard_soft_post$forecaster_loss)

losses <- round(losses, 3)

expect_true(length(unique(losses)) == length(losses))
# %%