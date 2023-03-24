library(splines2)
devtools::load_all()
devtools::load_all()

x <- 0:1000 / 1000
lambda <- 10
order <- 3
deg <- order - 1
n_inner <- 3 # Inner knots

knots <- make_knots2(n_inner, deg = deg)


spline_cpp <- make_basis_matrix2(x, knots, deg)

ts.plot(as.matrix(spline_cpp), col = 1:ncol(spline_cpp), lwd = 2)

p_spline_cpp <- splines2_periodic(
    x,
    knots,
    deg
)

ts.plot(as.matrix(p_spline_cpp), col = 1:ncol(p_spline_cpp), lwd = 2)
