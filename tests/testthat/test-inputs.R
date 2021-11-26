skip_if(debug_mode)

set.seed(1)

T <- 1000
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

# We expect inheritance from basis > smooth:
mod1 <- online(
    y = y, # Y as a vector
    tau = prob_grid,
    experts = experts,
    trace = FALSE
)

mod2 <- online(
    y = matrix(y), # Y as a vector
    tau = prob_grid,
    experts = experts,
    trace = FALSE
)

# Passting y as vec does not change results
expect_true(all(mod1$weights == mod2$weights))

