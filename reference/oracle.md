# Probabilistic Forecast Combination - Oracle

Returns predictions and weights calculated by numeric optimization. The
optimization is done in hindsight. This means all observations are used.

## Usage

``` r
oracle(y, experts, tau, affine = FALSE,
positive = FALSE, intercept = FALSE, debias = TRUE,
loss_function = "quantile", loss_parameter = 1, forget = 0)
```

## Arguments

- y:

  A numeric matrix of realizations. In probabilistic settings a matrix
  of dimension Tx1, in multivariate settings a TxD matrix. In the latter
  case, each slice of the expert's array gets evaluated using the
  corresponding column of the y matrix.

- experts:

  An array of predictions with dimension (Observations, Quantiles,
  Experts).

- tau:

  A numeric vector of probabilities.

- affine:

  Defines whether weights are summing to 1 or not. Defaults to FALSE.

- positive:

  Defines if a positivity constraint is applied to the weights. Defaults
  to FALSE.

- intercept:

  Determines if an intercept is added, defaults to FALSE. If true, a new
  first expert is added, always predicting 1.

- debias:

  Defines whether the intercepts weight is constrained or not. If TRUE
  (the default), the intercept weight is unconstrained. Only affects the
  results if affine and or positive is set to TRUE. If FALSE, the
  intercept is treated as an expert.

- loss_function:

  Either "quantile", "expectile" or "percentage".

- loss_parameter:

  Optional parameter scaling the power of the loss function.

- forget:

  Adds an exponential forgetting to the optimization. Past observations
  will get less influence on the optimization. Defaults to 0, which
  corresponds to no forgetting.

## Value

Returns weights and corresponding predictions. It is possible to
calculate the best convex combination of weights by setting affine and
positive to TRUE.

## Examples

``` r
if (FALSE) { # \dontrun{
T <- 50 # Observations
N <- 2 # Experts
P <- 9 # Quantiles
prob_grid <- 1:P / (P + 1)

y <- rnorm(n = T) # Realized
experts <- array(dim = c(T, P, N)) # Predictions
for (t in 1:T) {
    experts[t, , 1] <- qnorm(prob_grid, mean = -1, sd = 1)
    experts[t, , 2] <- qnorm(prob_grid, mean = 3, sd = sqrt(4))
}

model <- oracle(
    y = matrix(y),
    experts = experts
)
} # }
```
