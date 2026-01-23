# Predict method for online models

Calculates predictions based on new expert advice. This does not update
weights. If new observations are available use update instead. The
latter updates and weights and computes predictions.

## Usage

``` r
# S3 method for class 'online'
predict(object, new_experts, update_model = TRUE, ...)
```

## Arguments

- object:

  Object of class inheriting from 'online'

- new_experts:

  new expert predictions

- update_model:

  Defines whether the model object should be updated or not. If TRUE,
  new forecaster and expert predictions are appended onto the respective
  object items. Defaults to TRUE.

- ...:

  further arguments are ignored

## Value

`predict.online` produces an updated model object.
