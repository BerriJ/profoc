skip_if(debug_mode)

# %% Non periodic penalties with non equidistant knots
order <- 2
deg <- order - 1
J <- 5 # Number of inner knots. Total number of knots is J + 2*(order)
mu <- 0.4
sig <- 1

knots <- make_knots(J, deg = deg)

P1_cpp <- penalty(knots, order = deg + 1)[[1]]
D1 <- diff(diag(J + order), differences = 1)
P1_r <- t(D1) %*% D1

expect_equal(
    as.matrix(P1_cpp),
    P1_r
)

order <- 3
deg <- order - 1
knots <- make_knots(J, deg = deg)

P1_cpp <- penalty(knots, order = deg + 1)[[1]]
D1 <- diff(diag(J + order), differences = 1)
P1_r <- t(D1) %*% D1

expect_equal(
    as.matrix(P1_cpp),
    P1_r
)

P2_cpp <- penalty(knots, order = deg + 1)[[2]]
D2 <- diff(diag(J + order), differences = 2)
P2_r <- t(D2) %*% D2

expect_equal(
    as.matrix(P2_cpp),
    P2_r
)

order <- 4
deg <- order - 1
knots <- make_knots(J, deg = deg)

P2_cpp <- penalty(knots, order = deg + 1)[[3]]
D2 <- diff(diag(J + order), differences = 3)
P2_r <- t(D2) %*% D2

expect_equal(
    as.matrix(P2_cpp),
    P2_r
)
# %%

# %% Periodic penalties with equidistant knots
order <- 4
deg <- order - 1
mu <- 0.5
sig <- 1

n <- 6
Ip <- diag(n + 1)
Dp <- diff(Ip)[, -1]
Dp[1, n] <- -1

Dp2 <- t(Dp) %*% Dp
Dp3 <- t(Dp) %*% Dp2

PP1 <- t(Dp) %*% Dp
PP1
PP1_cpp <- penalty(knots, order, TRUE)[[1]]
expect_true(all(round(PP1_cpp, 12) == round(PP1, 12)))

PP2 <- t(Dp2) %*% Dp2
PP2_cpp <- penalty(knots, order, TRUE)[[2]]
expect_true(all(round(PP2_cpp, 12) == round(PP2, 12)))
PP3 <- t(Dp3) %*% Dp3
PP3_cpp <- penalty(knots, order, TRUE)[[3]]
expect_true(all(round(PP3_cpp, 12) == round(PP3, 12)))
# %%
