skip_if(debug_mode)
# %% Test lead for online
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

# Smoothing
boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    smooth_lambda = c(10, 1000),
    trace = FALSE
)

boa_smooth_post <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    smooth_lambda = c(10, 1000),
    smooth_ex_post = TRUE,
    trace = FALSE
)

expect_true(
    round(mean(boa_smooth$forecaster_loss), 6) !=
        round(mean(boa_smooth_post$forecaster_loss), 6)
)

boa_fs <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    fixed_share = 0.2,
    trace = FALSE
)

boa_fs_post <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    fixed_share = 0.2,
    fixed_share_ex_post = TRUE,
    trace = FALSE
)

expect_true(
    round(mean(boa_fs$forecaster_loss), 4) !=
        round(mean(boa_fs_post$forecaster_loss), 4)
)
# %%