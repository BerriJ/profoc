# Create a List of Hat Matrices

This function creates a list of hat matrices and the corresponding
parameters. It is used in
[`online()`](https://profoc.berrisch.biz/reference/online.md) to create
the hat matrices for penalized smoothing.

## Usage

``` r
make_hat_mats(
  x,
  n = length(x),
  mu = 0.5,
  sigma = 1,
  nonc = 0,
  tailw = 1,
  deg = 1,
  ndiff = 1.5,
  lambda = -Inf,
  periodic = FALSE,
  idx = NULL,
  params = NULL
)
```

## Arguments

- x:

  The predictor variable

- n:

  Number of knots

- mu:

  Beta distribution location parameter

- sigma:

  Beta distribution scale parameter

- nonc:

  Beta distribution noncentrality parameter

- tailw:

  Tailweight

- deg:

  Degree of splines

- ndiff:

  Sets the degree of the differencing matrix for creating the penalty

- lambda:

  Penalty parameter (higher values lead to higher penalty)

- periodic:

  Create periodic penalty

- idx:

  `make_hat_mats()` will create a grid containing all combinations of
  the parameters. If idx is set, this grid will be subsetted to the rows
  specified by idx.

- params:

  Instead of the arguments above, a grid (data.frame or named matrix) of
  parameters can be passed directly.
