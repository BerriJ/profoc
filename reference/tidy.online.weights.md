# Tidy the Weights of an Online object

`tidy` will transform the weights array of an online object into a
tibble that is better suited for plotting and analysis.

## Usage

``` r
# S3 method for class 'online.weights'
tidy(x, ...)
```

## Arguments

- x:

  The weights of an `online` object.

- ...:

  Not currently used.

## Value

A tibble with columns `t` `d` `p` `k` and `w` corresponding to the time,
marginals, probabilities, experts, and weights of the online-learning
computation.
