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
