# Simple setting of hat function

devtools::load_all()
devtools::load_all()

x <- 0:1000 / 100
lambda <- 10
order <- 2
deg <- order - 1
n_inner <- 3 # Inner knots
mu <- 0.3

knots <- c(0.00, 0.082, 0.23, 0.47, 1.00)





# install.packages("pbs")
# install.packages("splines2")

x <- 0:1000 / 1000
order <- 2
deg <- order - 1

knots <- c(0.00, 0.082, 0.23, 0.47, 1.00)

par(mfrow = c(1, 2))

ts.plot(
    pbs::pbs(x,
        # pbs add boundary knots automatically
        knots = knots[c(-1, -5)],
        degree = deg, intercept = TRUE
    ),
    col = seq_along(knots),
    lwd = 2
)

ts.plot(
    splines2::mSpline(x,
        knots = knots[c(-1, -5)],
        degree = deg,
        Boundary.knots = c(0, 1),
        periodic = TRUE,
        intercept = TRUE
    ),
    col = seq_along(knots),
    lwd = 2
)
