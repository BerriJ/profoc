Rcpp::compileAttributes()
devtools::document()
devtools::build()
devtools::check(cran = TRUE)
devtools::test()

# Release:

# Bump Version Number in Description
# Update the Date field in Description
# Update the Date in NEWS.md

# CRAN Release checklist:

devtools::spell_check()
devtools::check_rhub()
devtools::release()

# Merge submitted_to_cran into main, without squashing commits, wihtout extra merge commit
# Create a tag on main

# Potential library unload problems can be solved by specifying
# CXXFLAGS=-fno-gnu-unique in ~/.R/Makevars
# or
# PKG_CXXFLAGS=-fno-gnu-unique in src/Makevars

devtools::load_all()

devtools::build()
devtools::install()
library(profoc)
?online


# T <- 50 # Observations
# N <- 2 # Experts
# P <- 9 # Quantiles
# prob_grid <- 1:P / (P + 1)

# y <- rnorm(n = T) # Realized
# experts <- array(dim = c(T, P, N)) # Predictions
# for (t in 1:T) {
#     experts[t, , 1] <- qnorm(prob_grid, mean = -1, sd = 1)
#     experts[t, , 2] <- qnorm(prob_grid, mean = 3, sd = sqrt(4))
# }

# model <- online(
#     y = matrix(y),
#     experts = experts,
#     tau = prob_grid,
#     p_smooth_pr = list(lambda = 10)
# )

# print(model)
# plot(model)

# new_y <- matrix(rnorm(1)) # Realized
# new_experts <- experts[T, , , drop = FALSE]

# # Update will update the models weights etc if you provide new realizations
# model <- update(model, new_y = new_y, new_experts = new_experts)

# # Predict will expand \code{model$predictions} by default
# model <- predict(model, new_experts = new_experts, update_model = TRUE)