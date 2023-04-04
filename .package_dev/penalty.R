devtools::load_all()
devtools::load_all()

# %%
differences <- 1 # This script currently works for diff = 1 only
order <- 2
deg <- order - 1
J <- 3 # Number of inner knots. Total number of knots is J + 2*(order)
mu <- 0.5
sig <- 1
knots <- make_knots2(J, mu = mu, sig = sig, deg = deg)
# %%

# %% Non-periodic
# This is our target
D1D1_profoc <- as.matrix(penalty2(knots, order)[[1]])

# Manually:
h <- get_h(knots, order)
# D1 <- diff(diag(J + order), differences = differences)
# w_D1 <- diag(1 / h[, 1]) %*% D1
# D1D1_manual <- t(w_D1) %*% w_D1 # Penaltx matrix

# Using laplacian:
Adj <- matrix(0, nrow = J + order, ncol = J + order)
diag(Adj[-1, ]) <- 1
diag(Adj[, -1]) <- 1
plot(igraph::graph_from_adjacency_matrix(Adj, mode = "undirected"))
inc <- adjacency_to_incidence(Adj)
w_inc <- diag(1 / h[, 1]) %*% t(inc)
D1D1_inc <- -t(w_inc) %*% w_inc
diag(D1D1_inc) <- -diag(D1D1_inc)

D1D1_inc == D1D1_profoc
# %%

# %% Periodic and Equidistant

# For deg = 1 and first differences this is inner knots + 2 boundary knots
knots_sub <- knots[(differences + 1):(length(knots) - differences)]

h <- get_h(knots, order)[, 1]

Adj <- periodic_adjacency(J + order - 1)
inc_p <- adjacency_to_incidence(Adj)
w_inc_p <- diag(1 / h) %*% t(inc_p)
D1D1_inc_p <- -t(w_inc_p) %*% w_inc_p
diag(D1D1_inc_p) <- -diag(D1D1_inc_p)



plot(igraph::graph_from_adjacency_matrix(
    Adj,
    mode = "undirected"
))
image(D1D1_inc_p)

penalty_periodic(knots, order) == D1D1_inc_p * mean(h)^2
# %%


# %% Diff > 1
n <- 6
Ip <- diag(n + 1)
Dp <- diff(Ip)[, -(n + 1)]
Dp[n, 1] <- 1

Dp2 <- t(Dp) %*% Dp
Dp3 <- t(Dp2) %*% Dp2

Dp
Dp2
Dp3

knots <- make_knots2(5, mu = 0.5, deg = 1)
penalty(knots, order = 3)[[2]]
# %%
