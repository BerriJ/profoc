# Create a List of Basis Matrices

This function creates a list of basis matrices and the corresponding
parameters. It is used in
[`online()`](https://profoc.berrisch.biz/reference/online.md) to create
the basis matrices for basis smoothing.

## Usage

``` r
make_basis_mats(
  x,
  n = length(x),
  mu = 0.5,
  sigma = 1,
  nonc = 0,
  tailw = 1,
  deg = 1,
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

- periodic:

  Create periodic basis

- idx:

  `make_basis_mats()` will create a grid containing all combinations of
  the parameters. If idx is set, this grid will be subsetted to the rows
  specified by idx.

- params:

  Instead of the arguments above, a grid (data.frame or named matrix) of
  parameters can be passed directly.
