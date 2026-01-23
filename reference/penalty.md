# B-Spline penalty

This function calculates the B-Spline basis penalty. It follows the
procedure outlined in the paper by Zheyuan Li, Jiguo Cao, 2022 "General
P-Splines for Non-Uniform B-Splines"
[doi:10.48550/arXiv.2201.06808](https://doi.org/10.48550/arXiv.2201.06808)
. For equidistant knots it coincides with the usual penalty based on the
identitiy. For non-equidistant knots it is a weighted penalty with
respect to the knot distances. In addition to the above, we added the
possibility to calculate periodic penalties which are based on the
periodic differencing matrices.

## Usage

``` r
penalty(knots, order, periodic = FALSE, max_diff = 999L)
```

## Arguments

- knots:

  Vector of knots.

- order:

  Order of the Basis (degree + 1).

- periodic:

  Whether the penalties should be periodic or not.

- max_diff:

  Maximum difference order to calculate.

## Value

Returns a list of (order - 1) penalty matrices.

## Examples

``` r
if (FALSE) { # \dontrun{
# Equidistant knots with order 2
knots <- 1:10

P <- penalty(knots, order = 2)

print(P[[1]]) # First differences

# Non equidistant knots
knots <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)

P <- penalty(knots, order = 4)

print(P[[1]]) # First differences
print(P[[2]]) # Second differences
print(P[[3]]) # Third differences

# Periodic penalty for equidistant knots
oder <- 4
deg <- order - 1
knots <- 1:15

penalty(knots, order = order, periodic = TRUE)[[1]]
penalty(knots, order = order, periodic = TRUE)[[2]]
penalty(knots, order = order, periodic = TRUE)[[3]]
} # }
```
