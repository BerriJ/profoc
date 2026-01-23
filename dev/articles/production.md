# Production

## Using `online` in Production

This vignette explains the use of
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`update()`](https://rdrr.io/r/stats/update.html). These are the two
most important functions when using `profoc` in production. The
[`predict()`](https://rdrr.io/r/stats/predict.html) method is used to
combine new expert forecasts using the most recent combination weights.
This is useful if we combine new expert forecasts with the most recent
combination weights, but new observations have yet to be realized. At a
later point, [`update()`](https://rdrr.io/r/stats/update.html) can be
used to update the combination weights by evaluating the realized
observations. We assume that you followed the
[`vignette("profoc")`](https://profoc.berrisch.biz/dev/articles/profoc.md)
already. We will reuse the data and the model from there.

### Combining new expert predictions

First, we create new expert predictions:

``` r
new_experts <- experts[T, , , drop = FALSE]
```

The default behavior of
[`predict()`](https://rdrr.io/r/stats/predict.html) updates the
`combination` object. So, it can later be used to update the combination
weights as realized values emerge. That is,
[`predict()`](https://rdrr.io/r/stats/predict.html) expands
`combination$predictions` and returns the updated `combination`.

``` r
dim(combination$predictions)
#> [1] 32  1 99

# Predict will expand combination$predictions
combination <- predict(combination,
  new_experts = new_experts
)

dim(combination$predictions)
#> [1] 33  1 99
```

If you are only interested in the predictions, you can set
`update_model = FALSE`. In this case,
[`predict()`](https://rdrr.io/r/stats/predict.html) solely returns the
predictions:

``` r
predictions <- predict(combination,
  new_experts = new_experts,
  update_model = FALSE
)

dim(predictions)
#> [1]  1  1 99
```

### Updating the model weights

As new realizations emerge, we can update the combination weights. This
is done by [`update()`](https://rdrr.io/r/stats/update.html). That is,
[`update()`](https://rdrr.io/r/stats/update.html) expands
`combination$weights` and returns the updated `combination`.

``` r
# New observation
new_y <- matrix(rnorm(1))
```

``` r
dim(combination$weights)
#> [1] 33  1 99  2

# Model Update
combination <-
  update(combination,
    new_y = new_y
  )

dim(combination$weights)
#> [1] 34  1 99  2
```

### Summary on `predict()` and `update()`

As seen above, [`predict()`](https://rdrr.io/r/stats/predict.html) and
[`update()`](https://rdrr.io/r/stats/update.html) are closely related
and usually called sequentially. In an only setting, we want to
calculate the forecast (the combination) as soon as new expert
predictions emerge. For that, we can use
[`predict()`](https://rdrr.io/r/stats/predict.html). Later, as new
observations are realized, we can
[`update()`](https://rdrr.io/r/stats/update.html) the combination
weights.

We designed to also work in non-standard scenarios. So if, for example,
experts provide multi-step-ahead predictions, we can use
[`predict()`](https://rdrr.io/r/stats/predict.html) to combine all of
them using the most recent combination weights. Afterward, one or
multiple [`update()`](https://rdrr.io/r/stats/update.html) calls can be
used to update the combination weights as new observations are realized.
If we want to [`predict()`](https://rdrr.io/r/stats/predict.html) and
[`update()`](https://rdrr.io/r/stats/update.html) simultaneously, we can
do this. We can pass the new expert predictions and observations to
[`predict()`](https://rdrr.io/r/stats/predict.html). This will update
the combination weights and predictions with only one call to
[`predict()`](https://rdrr.io/r/stats/predict.html).
