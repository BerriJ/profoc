skip_if(debug_mode)
order <- 4
deg <- order - 1
a <- sqrt(2)

x <- 0:100 / 100

# Create desired knots:
knots <- make_knots(9, deg = deg)

# Create the b-spline basis using splines
old <- splines::splineDesign(knots, x, deg + 1, derivs = 0L, outer.ok = TRUE)

# Create the b-spline basis using splines2
new <- as.matrix(make_basis_matrix(x, knots, deg))
dimnames(new) <- NULL

expect_equal(old, new)

kstep <- 0.1
order <- 3
deg <- order - 1
a <- sqrt(2)
pstep <- 0.01

# Create desired knots:
knots <- make_knots(9, deg = deg)

# Create the b-spline basis using splines
old <- splines::splineDesign(knots, x, deg + 1, derivs = 0L, outer.ok = TRUE)

# Create the b-spline basis using splines2
new <- as.matrix(make_basis_matrix(x, knots, deg))
dimnames(new) <- NULL

expect_equal(old, new)

# Create knot sequence
order <- 2
deg <- order - 1
n_inner <- 3
mu <- 0.3
intercept <- TRUE

knots <- make_knots(n_inner, mu = mu, deg = deg)

K <- length(knots)

B <- splines2_basis(x, knots, deg, intercept = intercept)

# There should be K-order splines
expect_true(dim(B)[2] == K - order)

# Periodic Case
B <- splines2_basis(x, knots, deg, periodic = TRUE, intercept = intercept)
expect_true(dim(B)[2] == K - 2 * deg - 1)

order <- 4
deg <- order - 1

knots <- make_knots(n_inner, mu = mu, deg = deg)

K <- length(knots)

B <- splines2_basis(x, knots, deg, intercept = intercept)
expect_true(!all((B[1, ] - B[nrow(B), ]) == 0))

# There should be K-order splines
expect_true(dim(B)[2] == K - order)

# Periodic Case
B <- splines2_basis(x, knots, deg, periodic = TRUE, intercept = intercept)
expect_true(dim(B)[2] == K - 2 * deg - 1)

expect_true(all((B[1, ] - B[nrow(B), ]) == 0))
