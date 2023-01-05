skip_if(debug_mode)

x <- 1:9 / 10

deg <- 2

knots <- make_knots2(9, deg = deg)

D1 <- delta(knots, order = deg + 1)[[1]]

expect_equal(
    D1,
    diff(diag(12), differences = 1) * 10
)

D2 <- delta(knots, order = deg + 1)[[2]]

expect_equal(
    D2,
    diff(diag(12), differences = 2) * 100
)
