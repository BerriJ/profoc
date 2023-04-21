# remotes::install_github("igraph/rigraph")
# remotes::install_github("berrij/profoc")

devtools::load_all()
devtools::load_all()
library(igraph)


# %%
order <- 3
deg <- order - 1
n <- 3 # Number of inner knots. Total number of knots is n + 2*(order)
knots <- profoc:::make_knots(n, deg = deg)
# %%

# %% First differences non-periodic
## Valid for deg = 1, 2, 3 ...
D1 <- diff(diag(n + deg + 1), differences = 1)
D1D1 <- t(D1) %*% D1 # Penaltx matrix

# Adjacency matrix
Adj <- (D1D1 == -1) + 0

plot(graph_from_adjacency_matrix(Adj, mode = "undirected"))

# Example from SO: https://stackoverflow.com/a/22405325/9551847
# A <- matrix(
#     c(1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0),
#     nrow = 4,
#     ncol = 4
# )
# res <- adjacency_to_incidence(A)

inc <- adjacency_to_incidence(Adj)

penalty <- 2 * diag(apply(Adj, 1, sum)) - inc %*% t(inc)

penalty == D1D1

laplacian1 <- diag(n + deg + 1) * apply(Adj, 1, sum) - Adj
image(laplacian1)

D1D1 == laplacian1 # Laplacian is penalty matrix
as.matrix(profoc:::penalty(knots, order)[[1]]) == laplacian1 # Penalty above is the same as the one from profoc
# %%

# %% First differences periodic

## Valid for deg = 1, 2, 3 ...
B <- diag(n + deg + 1)
D1 <- diff(B, differences = 1)

D1D1 <- t(D1) %*% D1

## graphical way with adjacency matrix:
Adj <- (D1D1 == -1) + 0

Adj[nrow(Adj), 1] <- Adj[1, ncol(Adj)] <- 1

plot(graph_from_adjacency_matrix(Adj, mode = "undirected"))

laplacian1_p <- diag(n + deg + 1) * apply(Adj, 1, sum) - Adj
image(laplacian1_p)

inc <- adjacency_to_incidence(Adj)

penalty <- 2 * diag(apply(Adj, 1, sum)) - inc %*% t(inc)
penalty == laplacian1_p
# %%
