# Post Process Data from conline Class

This function works in conjunction with the conline class. After the
main learning task, it takes the output of the conline class and returns
an object suitable for, visualization, further, and deployment.
analysis.

## Usage

``` r
post_process_model(model_instance, names)
```

## Arguments

- model_instance:

  An instance of conline.

- names:

  A named list with dimnames of `y` and `experts`.
