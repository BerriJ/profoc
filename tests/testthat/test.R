# Experts
N <- 2
# Observations
T <- 100
# Size of probability grid
P <- 99
prob_grid <- 1:P / (P + 1)

# Realized observations
y <- rnorm(T)

# Deviation of the experts
dev <- c(-1, 3)
experts_sd <- c(1, sqrt(4))

# Expert predictions
experts <- array(dim = c(T, P, N))

for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = dev[1], sd = experts_sd[1])
    experts[t, , 2] <- qnorm(prob_grid, mean = dev[2], sd = experts_sd[2])
}

expect_type(profoc(
    y = matrix(y),
    experts = experts
), "list")