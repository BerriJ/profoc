# remotes::install_github("igraph/rigraph")
# remotes::install_github("berrij/profoc")

devtools::load_all()
devtools::load_all()

# %%
order <- 3
deg <- order - 1
n <- 3 # Number of inner knots. Total number of knots is n + 2*(order)
knots <- profoc:::make_knots2(n, deg = deg)
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

penalty <- 2 * diag(apply(Adj, 1, sum)) - t(inc) %*% inc

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

penalty <- 2 * diag(apply(Adj, 1, sum)) - t(inc) %*% inc
penalty == laplacian1_p

# How to get D1 from the laplacian matrix? We need to adjust it using: solve(W) %*% D1 to calculate the penalty matrix for non-equidistant knots
# See eq. (18) https://arxiv.org/pdf/2303.10019v1.pdf
# %%

## %% Second differences non-periodic
D2 <- diff(diag(n + deg + 1), differences = 2)
D2D2 <- t(D2) %*% D2

# This one was created by "reverse engineering" the penalty matrix
Adj <- matrix(c(
    0, 2, -1, 0, 0, 0,
    2, 0, 4, -1, 0, 0,
    -1, 4, 0, 4, -1, 0,
    0, -1, 4, 0, 4, -1,
    0, 0, -1, 4, 0, 2,
    0, 0, 0, -1, 2, 0
), byrow = TRUE, 6, 6)

# Smarter way would be to use the penalty matrix and L = D - A
# 1. Take the penalty matrix
Adj_ <- D2D2
# 2. Get the negative adjacency by seting the diagonal to 0
diag(Adj_) <- 0
# 3. Switch signs
Adj_ <- -Adj_

Adj_ == Adj

# Why the minus sign here? Without it, the graph makes no sense?
plot(graph_from_adjacency_matrix(-Adj_, mode = "undirected"))

laplacian2 <- diag(n + deg + 1) * apply(Adj, 1, sum) - Adj
image(laplacian2)

D2D2 == laplacian2 # Laplacian is penalty matrix
as.matrix(profoc:::penalty(knots, order)[[2]]) == laplacian2 # Penalty above is the same as the one from profoc
# %%

## %% Second differences periodic
D2 <- diff(diag(n + deg + 1), differences = 2)
D2D2 <- t(D2) %*% D2

Adj <- matrix(c(
    0, 4, -1, 0, -1, 4,
    4, 0, 4, -1, 0, -1,
    -1, 4, 0, 4, -1, 0,
    0, -1, 4, 0, 4, -1,
    -1, 0, -1, 4, 0, 4,
    4, -1, 0, -1, 4, 0
), byrow = TRUE, 6, 6)

plot(graph_from_adjacency_matrix(-Adj, mode = "directed"))

laplacian2_p <- diag(n + deg + 1) * apply(Adj, 1, sum) - Adj
image(laplacian2_p)
# %%
