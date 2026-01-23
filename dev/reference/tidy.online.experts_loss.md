# Tidy the Experts' losses of an Online object

`tidy` will transform the `experts_loss` array of an online object into
a tibble that is better suited for plotting and analysis.

## Usage

``` r
# S3 method for class 'online.experts_loss'
tidy(x, ...)
```

## Arguments

- x:

  The experts_loss of an `online` object.

- ...:

  Not currently used.

## Value

A tibble with columns `t` `d` `p` `k` and `w` corresponding to the time,
marginals, probabilities, and experts_loss of the online-learning
computation.
