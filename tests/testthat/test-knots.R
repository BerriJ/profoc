skip_if(debug_mode)

# Defaults should create equidistant grid and identity
x <- 1:9 / 10
deg <- 1

knots <- make_knots(n = length(x))

expect_true(all(round(diff(knots), 10) == 0.1))

B <- make_basis_mats(x)$basis[[1]]

expect_true(all(round(as.matrix(B), 10) == diag(rep.int(1, length(x)))))

# n-1 should create a vector of ones
knots <- make_knots(n = -1)
expect_null(knots)

expect_true(all(matrix(rep.int(1, length(x))) ==
    as.matrix(make_basis_mats(x, n = -1)$basis[[1]])))

# Now check different specifications and check expected dims
deg <- 1
n <- 4

B1 <- make_basis_mats(x, n = n, deg = 1)$basis[[1]]
B2 <- make_basis_mats(x, n = n, deg = 2)$basis[[1]]
B3 <- make_basis_mats(x, n = n, deg = 3)$basis[[1]]

expect_true(sum(as.matrix(B1) != as.matrix(B2)) == 26)
expect_true(sum(as.matrix(B2) != as.matrix(B3)) == 36)
expect_true(sum(as.matrix(B1) != as.matrix(B3)) == 36)

# Check if dist parameters change matrix
B1 <- make_basis_mats(x, n = n, mu = 0.5, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, mu = 0.6, deg = 3)$basis[[1]]
expect_true(any(as.matrix(B1) != as.matrix(B2)))

B1 <- make_basis_mats(x, n = n, sigma = 1, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, sigma = 1.2, deg = 3)$basis[[1]]
expect_true(any(as.matrix(B1) != as.matrix(B2)))

B1 <- make_basis_mats(x, n = n, nonc = 1, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, nonc = 1.2, deg = 3)$basis[[1]]
expect_true(any(as.matrix(B1) != as.matrix(B2)))

B1 <- make_basis_mats(x, n = n, nonc = 1, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, nonc = 1.2, deg = 3)$basis[[1]]
expect_true(any(as.matrix(B1) != as.matrix(B2)))

B1 <- make_basis_mats(x, n = n, tailw = 1, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, tailw = 1.2, deg = 3)$basis[[1]]
expect_true(any(as.matrix(B1) != as.matrix(B2)))

# Control
B1 <- make_basis_mats(x, n = n, tailw = 1, deg = 3)$basis[[1]]
B2 <- make_basis_mats(x, n = n, tailw = 1, deg = 3)$basis[[1]]
expect_false(any(as.matrix(B1) != as.matrix(B2)))


# Check if dist parameters change hat matrix
H1 <- make_hat_mats(x, n = n, lambda = 5, mu = 0.5, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, mu = 0.6, deg = 3)$hat[[1]]
expect_true(any(as.matrix(H1) != as.matrix(H2)))

H1 <- make_hat_mats(x, n = n, lambda = 5, sigma = 1, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, sigma = 1.2, deg = 3)$hat[[1]]
expect_true(any(as.matrix(H1) != as.matrix(H2)))

H1 <- make_hat_mats(x, n = n, lambda = 5, nonc = 1, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, nonc = 1.2, deg = 3)$hat[[1]]
expect_true(any(as.matrix(H1) != as.matrix(H2)))

H1 <- make_hat_mats(x, n = n, lambda = 5, nonc = 1, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, nonc = 1.2, deg = 3)$hat[[1]]
expect_true(any(as.matrix(H1) != as.matrix(H2)))

H1 <- make_hat_mats(x, n = n, lambda = 5, tailw = 1, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, tailw = 1.2, deg = 3)$hat[[1]]
expect_true(any(as.matrix(H1) != as.matrix(H2)))

# Control
H1 <- make_hat_mats(x, n = n, lambda = 5, tailw = 1, deg = 3)$hat[[1]]
H2 <- make_hat_mats(x, n = n, lambda = 5, tailw = 1, deg = 3)$hat[[1]]
expect_false(any(as.matrix(H1) != as.matrix(H2)))

# Check online estimation:

# Experts
N <- 2
# Observations
T <- 1000
D <- 2
# Size of probability grid
P <- 99
prob_grid <- 1:P / (P + 1)

# Realized observations
y <- matrix(rnorm(T), nrow = T, ncol = 2)

dev <- c(-2, 2)
experts_sd <- c(1, 2)

# Expert predictions
experts <- array(dim = c(T, D, P, N))

for (t in 1:T) {
    experts[t, , , 1] <- qnorm(prob_grid, mean = dev[1], sd = experts_sd[1])
    experts[t, , , 2] <- qnorm(prob_grid, mean = dev[2], sd = experts_sd[2])
}

foo <- online(
    y = y,
    experts = experts,
    tau = prob_grid,
    trace = FALSE
)

# plot(foo$weights[T, 1, , 1], type = "l")

# Test b_smooth_pr
foo2 <- online(
    y = y,
    experts = experts,
    tau = prob_grid,
    b_smooth_pr = list(
        knots = 20,
        deg = 3
    ),
    trace = FALSE
)

# lines(foo2$weights[T, 1, , 1], type = "l", col = 2)

# Test p_smooth_pr
foo3 <- online(
    y = y,
    experts = experts,
    tau = prob_grid,
    p_smooth_pr = list(
        knots = 20,
        mu = 0.4,
        deg = 3,
        bdiff = 1.3,
        lambda = 5
    ),
    trace = FALSE
)

# lines(foo3$weights[T, 1, , 1], type = "l", col = 2)

# Test b_smooth_mv note that we use more knots > D which is generally possible
foo2 <- online(
    y = y,
    experts = experts,
    tau = prob_grid,
    b_smooth_mv = list(
        knots = 20,
        mu = 0.1,
        sigma = 5,
        nonc = 2,
        tailweight = 3,
        deg = 1
    ),
    trace = FALSE
)

# Test p_smooth_mv
foo2 <- online(
    y = y,
    experts = experts,
    tau = prob_grid,
    p_smooth_mv = list(
        knots = 20,
        mu = 0.1,
        sigma = 5,
        nonc = 2,
        tailweight = 3,
        deg = 1,
        bdiff = 1.3,
        lambda = 5
    ),
    trace = FALSE
)
