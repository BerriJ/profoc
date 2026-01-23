# Create experts list to be used in conline class

This function works in conjunction with the conline class. It takes a
matrix of experts and a matrix of outcomes and returns a list of experts
which fulfills all properties that are needed for passing it to the an
instance of conline.

## Usage

``` r
init_experts_list(experts, y, output_with_names = FALSE)
```

## Arguments

- experts:

  array of predictions with dimension T x D x P x K (Observations x
  Variables x Quantiles x Experts) or T x D x K or T x P x K.

- y:

  A matrix of outcomes with dimension T x D.

- output_with_names:

  Defaults to FALSE. If TRUE, the function returns a list with the
  experts list, the names of the variables (dnames) and the names of the
  experts (enames).
