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
    experts = experts,
    tau = prob_grid,
    trace = FALSE
)
# %%

# %% Estimate the first partial model
model_partial_1 <- online(
    y = matrix(y)[1:10],
    experts = experts,
    tau = prob_grid,
    trace = FALSE
)

model_partial_1_preds <- model_partial_1$predictions[11:20, , ]
# %%

# %% Update using same number of y and experts
model_partial_2 <- online(
    y = matrix(y[1:10]),
    tau = prob_grid,
    experts = experts[1:10, , , drop = FALSE],
    trace = FALSE
)
model_partial_2 <- predict(
    model_partial_2,
    new_experts = experts[11:20, , , drop = FALSE]
)

expect_true(all(model_partial_1$predictions == model_partial_2$predictions))
identical(model_partial_2, model_full)
# %%

# %% Asymetric updating
model_full_updated <- update(
    model_partial_1,
    new_y = matrix(y[11:20])
)

identical(model_full_updated, model_full)
# %%

# %% Don't update model
model_partial_3 <- online(
    y = matrix(y)[1:10],
    experts = experts[1:10, , ],
    tau = prob_grid,
    trace = FALSE
)
model_partial_3_copy <- model_partial_3
preds <- predict.online(model_partial_3, experts[11:20, , , drop = FALSE], update_model = FALSE)
expect_true(identical(model_partial_3, model_partial_3_copy))

expect_true(all(preds[1:10, , ] == model_partial_1_preds))

# %%
