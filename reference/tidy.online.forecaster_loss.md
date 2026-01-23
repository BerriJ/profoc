# Tidy the Experts' losses of an Online object

`tidy` will transform the \`forecaster_loss“ array of an online object
into a tibble that is better suited for plotting and analysis.

## Usage

``` r
# S3 method for class 'online.forecaster_loss'
tidy(x, ...)
```

## Arguments

- x:

  The forecaster_loss of an `online` object.

- ...:

  Not currently used.

## Value

A tibble with columns `t` `d` `p` `k` and `w` corresponding to the time,
marginals, probabilities, and forecaster_loss of the online-learning
computation.
