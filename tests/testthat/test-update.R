skip_if(debug_mode)

# %% Simulate data
set.seed(1)

T <- 20
N <- 2
P <- 9
prob_grid <- 1:P / (P + 1)

mean_y <- 0
sd_y <- 5

# Realized observations
y <- rnorm(n = T)

# Expert predictions
experts <- array(dim = c(T, P, N))
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = -5, sd = 2)
    experts[t, , 2] <- qnorm(prob_grid, mean = 5, sd = 2)
}

# %%

# %% Estimate the full model
model_full <- online(
    y = matrix(y),
    tau = prob_grid,
    experts = experts,
    trace = FALSE
)
# %%

# %% Update using same number of y and experts
model_partial <- online(
    y = matrix(y[1:10]),
    tau = prob_grid,
    experts = experts[1:10, , , drop = FALSE],
    trace = FALSE
)
model_partial <- update(
    object = model_partial,
    new_y = matrix(y[11:15]), new_experts = experts[11:15, , , drop = FALSE]
)
# %%

# %%
model_partial <- update(model_partial, new_y = matrix(y[16:20]), new_experts = experts[16:20, , , drop = FALSE])
# %%

# %% Models should now be identical
expect_true(identical(model_partial, model_full))
# %%

# %% Asymetric updating
model_partial <- online(
    y = matrix(y[1:5]),
    experts = experts,
    tau = prob_grid,
    trace = FALSE
)

model_partial <- update(
    model_partial,
    new_y = matrix(y[6:10])
)

model_partial <- update(
    model_partial,
    new_y = matrix(y[11:20])
)
# %%

# %% Models should now be identical
expect_true(identical(model_partial, model_full))
# %%

# %% Throw an error when insufficient amount of expert predictions is supplied
expect_error(
    update(model_partial,
        new_y = matrix(y[11:20]),
    ), "Number of provided expert predictions has to match or exceed observations."
)
# %%