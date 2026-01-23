# Tidy the Predictions of an Online object

`tidy` will transform the `predictions` array of an online object into a
tibble that is better suited for plotting and analysis.

## Usage

``` r
# S3 method for class 'online.predictions'
tidy(x, ...)
```

## Arguments

- x:

  The predictions of an `online` object.

- ...:

  Not currently used.

## Value

A tibble with columns `t` `d` `p` `k` and `w` corresponding to the time,
marginals, probabilities, and predictions of the online-learning
computation.
