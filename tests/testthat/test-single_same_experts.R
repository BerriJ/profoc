skip_if(debug_mode)
# %% Test online "combination" of a single expert
set.seed(1)

mod <- online(
    y = array(rnorm(30),
        dim = c(5, 3)
    ), array(rnorm(30),
        dim = c(5, 3, 1)
    ),
    tau = .5,
    trace = FALSE
)

expect_true(all(mod$weights == 1))
# %%


# %% Test online "combination" of two experts that are the same
set.seed(1)

experts <- array(NA,
    dim = c(5, 3, 2)
)

experts[, , 1] <- array(rnorm(30),
    dim = c(5, 3)
)
experts[, , 2] <- experts[, , 1]

mod <- online(
    y = array(rnorm(30),
        dim = c(5, 3)
    ), experts,
    tau = .5,
    trace = FALSE
)

expect_true(all(mod$weights == 0.5))
# %%
