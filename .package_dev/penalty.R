devtools::load_all()
devtools::load_all()

# %%
differences <- 3
order <- 4
deg <- order - 1
J <- 5 # Number of inner knots. Total number of knots is J + 2*(order)
mu <- 0.5
sig <- 1
knots <- make_knots2(J, mu = mu, sig = sig, deg = deg)

image(as.matrix(penalty(knots, order)[[differences]]))
image(as.matrix(penalty(knots, order, periodic = TRUE)[[differences]]))
# %%

# %% Diff > 1
devtools::load_all()

# Non Periodic
D1 <- diff(diag(J + deg + 1), differences = 1)
D1D1 <- t(D1) %*% D1 # Penaltx matrix

# Periodic
n <- J + 1
Ip <- diag(n + 1)
Dp <- diff(Ip)[, -1]
Dp[1, n] <- -1

Dp2 <- t(Dp) %*% Dp
Dp3 <- t(Dp) %*% Dp2

PP1 <- t(Dp) %*% Dp
PP1
PP1_cpp <- penalty(knots, order, TRUE)[[1]]
round(PP1_cpp, 3) == round(PP1, 3)

PP2 <- t(Dp2) %*% Dp2
PP2_cpp <- penalty(knots, order, TRUE)[[2]]
round(PP2_cpp, 3) == round(PP2, 3)
PP3 <- t(Dp3) %*% Dp3
PP3_cpp <- penalty(knots, order, TRUE)[[3]]
round(PP3_cpp, 3) == round(PP3, 3)
# %%

# %% Weighting the penalty w.r.t. knots

# This is the same as for weighting the periodic m-splines
idx_inner <- (order + 1):(length(knots) - order)
idx_bounds <- c(order, length(knots) - order + 1)
bounds <- knots[c(order, length(knots) - order + 1)]

knots_seq <- c(
    knots[idx_bounds[1]:idx_bounds[2]], # Inner and bounds
    knots[head(idx_inner, deg)] + knots[idx_bounds[2]] - knots[idx_bounds[1]]
)

K_ <- length(knots_seq)

n <- K_ - 2

diff(knots_seq[1:(K_ - 1)],
    lag = order - 1,
    differences = 1
)
knots_seq[(deg + 1):(K_ - 1)] - knots_seq[1:(K_ - order)]
# %%

# %%
par(mfrow = c(3, 5))
for (i in 2:4) {
    for (m in c(0.2, 0.4, 0.5, 0.6, 0.8)) {
        order <- i
        deg <- order - 1
        mu <- m
        n <- 10
        knots <- make_knots2(n, mu = mu, deg = deg)

        P <- as.matrix(penalty(knots, order, periodic = TRUE)[[order - 1]])

        # P <- as.matrix(penalty(knots, order)[[1]])

        PPPP <- matrix(NA, nrow = 2 * ncol(P), ncol = 2 * ncol(P))
        PPPP[1:ncol(P), 1:ncol(P)] <- P
        PPPP[1:ncol(P) + ncol(P), 1:ncol(P)] <- P
        PPPP[1:ncol(P), 1:ncol(P) + ncol(P)] <- P
        PPPP[1:ncol(P) + ncol(P), 1:ncol(P) + ncol(P)] <- P
        image(PPPP, main = paste("Deg=", deg, ", mu=", mu))
    }
}
# %%
