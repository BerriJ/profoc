# Add code to debug below
# Set desired breakpoints in .cpp files
# Start debugger with "Debug RCPP" in VS-Code

test_that("Debug", {
    skip_if_not(debug_mode)

    # Probabilistic example
    set.seed(1)

    # Quantiles ####################################################################

    # Experts
    N <- 2
    # Observations
    T <- 50
    # Size of probability grid
    P <- 9
    prob_grid <- 1:P / (P + 1)

    # Realized observations
    y <- rnorm(T)

    # Deviation of the experts
    # dev <- c(-1, 3)
    # experts_sd <- c(1, sqrt(4))

    dev <- c(-1, 3)
    experts_sd <- c(1, 4)

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
        lambda = 5,
        ndiff = 2,
        deg = 2
    )
})