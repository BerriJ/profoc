# Update method for online models

Continues learning using new observations and new expert advice.

## Usage

``` r
# S3 method for class 'online'
update(object, new_y, new_experts = NULL, trace = FALSE, ...)
```

## Arguments

- object:

  Object of class inheriting from 'online'

- new_y:

  new observations

- new_experts:

  new expert predictions. This must be left unspecified

- trace:

  If a progress bar shall be shown. Defaults to FALSE if the model
  already contains the expert predictions corresponding to new_y.

- ...:

  further arguments are ignored

## Value

`update.online` produces an updated model object.
