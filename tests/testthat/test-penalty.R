skip_if(debug_mode)

deg <- 1
knots <- make_knots2(9, deg = deg)

P1_cpp <- penalty(knots, order = deg + 1)[[1]]
D1 <- diff(diag(11), differences = 1)
P1_r <- t(D1) %*% D1

expect_equal(
    as.matrix(P1_cpp),
    P1_r
)

deg <- 2
knots <- make_knots2(9, deg = deg)

P1_cpp <- penalty(knots, order = deg + 1)[[1]]
D1 <- diff(diag(12), differences = 1)
P1_r <- t(D1) %*% D1

expect_equal(
    as.matrix(P1_cpp),
    P1_r
)

P2_cpp <- penalty(knots, order = deg + 1)[[2]]
D2 <- diff(diag(12), differences = 2)
P2_r <- t(D2) %*% D2

expect_equal(
    as.matrix(P2_cpp),
    P2_r
)

deg <- 3
knots <- make_knots2(9, deg = deg)

P2_cpp <- penalty(knots, order = deg + 1)[[3]]
D2 <- diff(diag(13), differences = 3)
P2_r <- t(D2) %*% D2

expect_equal(
    as.matrix(P2_cpp),
    P2_r
)

n <- 6
Ip <- diag(n + 1)
Dp <- diff(Ip)[, -1]
Dp[1, n] <- -1

Dp2 <- t(Dp) %*% Dp
Dp3 <- t(Dp) %*% Dp2

PP1 <- t(Dp) %*% Dp
PP1
PP1_cpp <- penalty_periodic2(knots, order)[[1]]
round(PP1_cpp, 3) == round(PP1, 3)

PP2 <- t(Dp2) %*% Dp2
PP2_cpp <- penalty_periodic2(knots, order)[[2]]
round(PP2_cpp, 3) == round(PP2, 3)
PP3 <- t(Dp3) %*% Dp3
PP3_cpp <- penalty_periodic2(knots, order)[[3]]
round(PP3_cpp, 3) == round(PP3, 3)
