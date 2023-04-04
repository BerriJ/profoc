devtools::load_all()
devtools::load_all()

par(mfrow = c(3, 5))
for (i in 2:4) {
    for (m in c(0.2, 0.4, 0.5, 0.6, 0.8)) {
        order <- i
        deg <- order - 1
        mu <- m
        n <- 10
        knots <- make_knots2(n, mu = mu, deg = deg)

        P <- penalty_periodic(knots, order)

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
order <- 2
deg <- order - 1
n <- 3 # Number of inner knots. Total number of knots is n + 2*(order)
mu <- 0.4
knots <- profoc:::make_knots2(n, mu = mu, deg = deg)

D1 <- diff(diag(n + deg + 1), differences = 1)
D1D1 <- t(D1) %*% D1 # Penaltx matrix

K <- length(knots)

# TODO: hardcore this 2 in cpp
adj <- periodic_adjacency(
    K - 2 * order + 1
)
inc <- adjacency_to_incidence(adj)

diff(knots[2:(K - 1)], lag = order - 1, differences = 1)
knots[(order + 1):(K - 1)] - knots[2:(K - order)]
penalty_periodic(knots, order)

penalty_periodic2(knots, order)

image(penalty_periodic2(knots, order))
# %%

# %% This is the same as for weighting the periodic m-splines
idx_inner <- (order + 1):(length(knots) - order)
idx_bounds <- c(order, length(knots) - order + 1)
bounds <- knots[c(order, length(knots) - order + 1)]

knots_seq <- c(
    knots[idx_bounds[1]:idx_bounds[2]], # Inner and bounds
    knots[head(idx_inner, deg)] + knots[idx_bounds[2]]
)

K_ <- length(knots_seq)

n <- K_ - 2

diff(knots_seq[1:(K_ - 1)],
    lag = order - 1,
    differences = 1
)
knots_seq[(deg + 1):(K_ - 1)] - knots_seq[1:(K_ - order)]
# %%
