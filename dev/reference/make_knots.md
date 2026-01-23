# Create a vector of knots for splines

This function creates a knot vector for splines. The knots are
distributed according to a beta distribution. The first input defines
the number of inner knots. The total number of knots is `n + 2 * order`.

## Usage

``` r
make_knots(n, mu = 0.5, sig = 1, nonc = 0, tailw = 1, deg = 1)
```

## Arguments

- n:

  Number of knots

- mu:

  Beta distribution location parameter

- sig:

  Beta distribution scale parameter

- nonc:

  Beta distribution noncentrality parameter

- tailw:

  Tailweight

- deg:

  Degree of splines
