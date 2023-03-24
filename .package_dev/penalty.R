devtools::load_all()
devtools::load_all()

# %%
differences <- 1

order <- 3
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
D1 <- diff(diag(J + order), differences = differences)
w_D1 <- diag(1 / h[, 1]) %*% D1
D1D1_manual <- t(w_D1) %*% w_D1 # Penaltx matrix

# Using laplacian:
Adj <- matrix(0, nrow = J + order, ncol = J + order)
diag(Adj[-1, ]) <- 1
diag(Adj[, -1]) <- 1
plot(igraph::graph_from_adjacency_matrix(Adj, mode = "undirected"))
inc <- adjacency_to_incidence(Adj)
w_inc <- diag(1 / h[, 1]) %*% inc
D1D1_inc <- -t(w_inc) %*% w_inc
diag(D1D1_inc) <- -diag(D1D1_inc)

D1D1_inc == D1D1_profoc
# %%

# %% Periodic

# Using laplacian:
h <- diff_cpp2(knots[(differences + 1):(length(knots) - differences)], order - differences, 1) / (order - differences)

h <- get_h(knots, order)

total <- length(knots)
outer <- 2 * order
inner_knots <- knots[(outer / 2):(total - outer / 2 + 1)]

Adj <- matrix(0, nrow = J + order, ncol = J + order)
diag(Adj[-1, ]) <- 1
diag(Adj[, -1]) <- 1
Adj[nrow(Adj), 1] <- 1
Adj[1, ncol(Adj)] <- 1
plot(igraph::graph_from_adjacency_matrix(Adj, mode = "undirected"))
inc <- adjacency_to_incidence(Adj)
w_inc <- diag(1 / h[, 1]) %*% inc
D1D1_inc <- -t(w_inc) %*% w_inc
diag(D1D1_inc) <- -diag(D1D1_inc)

D1D1_inc == D1D1_profoc
# %%
