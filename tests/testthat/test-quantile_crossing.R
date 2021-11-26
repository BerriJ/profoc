skip_if(debug_mode)

set.seed(1)

T <- 1000
N <- 2
P <- 99
prob_grid <- 1:P / (P + 1)

# Deviation of the experts
dev <- c(-1, 3)
experts_sd <- c(1, sqrt(4))

E1 <- qnorm(prob_grid, mean = dev[1], sd = experts_sd[1])
E2 <- qnorm(prob_grid, mean = dev[2], sd = experts_sd[2])


# Realized observations
y <- rnorm(T)

# Expert predictions
experts <- array(dim = c(T, P, N))

for (t_ in 1:T) {
    experts[t_, , 1] <- qnorm(prob_grid,
        mean = dev[1],
        sd = experts_sd[1]
    )
    experts[t_, , 2] <- qnorm(prob_grid,
        mean = dev[2],
        sd = experts_sd[2]
    )
}

boa_smooth_crossed <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    allow_quantile_crossing = TRUE,
    trace = FALSE
)

boa_smooth <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    allow_quantile_crossing = FALSE,
    trace = FALSE
)

# We expect that quantile crossing happens and therefore losses are different:
expect_true(
    sum(boa_smooth_crossed$forecaster_loss) != sum(boa_smooth$forecaster_loss)
)

# We expect that the "sorted loss" is less than the "unsorted loss"
expect_lt(
    sum(boa_smooth$forecaster_loss), sum(boa_smooth_crossed$forecaster_loss)
)

# Multivariate setting (sorting will be disabled automatically)
# T <- 1000
# N <- 2
# P <- 2
# prob_grid <- c(0.5, 0.5)

# mean_y <- 0
# sd_y <- 5

# # Realized observations
# y <- matrix(c(rnorm(n = T), rnorm(n = T)), ncol = 2)

# # Expert predictions
# experts <- array(dim = c(T, P, N))
# for (t in 1:T) {
#     experts[t, , 1] <- qnorm(prob_grid, mean = -5, sd = 2)
#     experts[t, , 2] <- qnorm(prob_grid, mean = 5, sd = 2)
# }

# boa_smooth_crossed <- online(
#     y = y,
#     tau = prob_grid,
#     experts = experts,
#     allow_quantile_crossing = TRUE,
#     trace = FALSE
# )

# expect_warning(
#     boa_smooth <- online(
#         y = y,
#         tau = prob_grid,
#         experts = experts,
#         allow_quantile_crossing = FALSE,
#         trace = FALSE
#     )
# )


# # We expect that the "sorted loss" is less than the "unsorted loss"
# expect_equal(
#     sum(boa_smooth$forecaster_loss), sum(boa_smooth_crossed$forecaster_loss)
# )