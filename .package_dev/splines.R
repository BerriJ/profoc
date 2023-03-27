library(splines2)
devtools::load_all()
devtools::load_all()

x <- 0:1000 / 1000
lambda <- 10
order <- 2
deg <- order - 1
n_inner <- 3 # Inner knots

knots <- make_knots2(n_inner, deg = deg)

spline_cpp <- make_basis_matrix2(
    x,
    knots,
    deg
)

ts.plot(
    as.matrix(spline_cpp),
    col = 1:ncol(spline_cpp), lwd = 2
)

spline_cpp_p <- splines2_periodic(
    x,
    knots,
    deg
)

ts.plot(
    as.matrix(spline_cpp),
    col = 1:ncol(spline_cpp), lwd = 2
)










grid <- expand.grid(
    periodic = c(FALSE, TRUE),
    n_inner = c(3, 8),
    mu = c(0.5, 0.3)
)

par(mfrow = c(4, 2))

for (i in 1:nrow(grid)) {
    knots <- make_knots2(n = grid[i, "n_inner"], mu = grid[i, "mu"], deg = deg)
    if (!grid[i, "periodic"]) {
        spline_cpp <- make_basis_matrix2(x, knots, deg)
    } else {
        spline_cpp <- splines2_periodic(
            x,
            knots,
            deg
        )
    }
    ts.plot(
        as.matrix(spline_cpp),
        col = 1:ncol(spline_cpp), lwd = 2,
        main = paste0("mu = ", grid[i, "mu"], ", inner = ", grid[i, "n_inner"], ", periodic = ", grid[i, "periodic"])
    )
}

spline_cpp <- make_basis_matrix2(x, knots, deg)

ts.plot(
    as.matrix(spline_cpp),
    col = 1:ncol(spline_cpp), lwd = 2,
    main = "First order B-splines"
)


p_spline_cpp <- splines2_periodic(
    x,
    knots,
    deg
)

ts.plot(as.matrix(p_spline_cpp), col = 1:ncol(p_spline_cpp), lwd = 2)
