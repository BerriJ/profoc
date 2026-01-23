# Using the C++ Interface

## Introduction

All major parts of
[`online()`](https://profoc.berrisch.biz/reference/online.md) are
implemented in C++ for speed. Usually, this comes at the cost of
flexibility. However, the profoc package exposes a C++ class `conline`
that allows you to gain fine grained control over objects.  
[`online()`](https://profoc.berrisch.biz/reference/online.md) wraps this
class and provides a convenient interface for the most common use cases.
However, if you need to alter object initialization (i.e. provide custom
basis / hat matrices for smoothing) you can use the C++ class directly
from R. This vignette shows how to do this.

Note that we will reuse the data from
[`vignette("profoc")`](https://profoc.berrisch.biz/articles/profoc.md).

## Online learning with `conline`

First, we need to create a new instance of the c++ class. This can be
done by calling `new(conline)`.

``` r
library(profoc)
model <- new(conline)
```

Now we need to pass the data to the class instance. The whole list of
accessible field can be printed with `names(model)`. Most of them have
defaults.

``` r
model$y <- y
tau <- 1:P / (P + 1)
model$tau <- tau
```

The experts array is a bit more complicated. C++ expects us to pass a
list of arrays. Thereby, the list itself must have dimension `Tx1` and
the elements of the list (the arrays) `D x P x K`. For convenience we
can use
[`init_experts_list()`](https://profoc.berrisch.biz/reference/init_experts_list.md)
to create such a list from our experts array. Note that we must pass the
true observations as well. They are used to detect whether the data is
univariate (`T x 1` matrix) or multivariate (`T x D` matrix).

``` r
experts_list <- init_experts_list(experts, y)
model$experts <- experts_list
```

Now suppose we want to alter the smoothing behavior across quantiles. We
start by creating a new hat matrix.

``` r
hat <- make_hat_mats(
  x = tau,
  mu = 0.2, # Put more knots in the lower tail
  periodic = TRUE
)
str(hat)
#> List of 2
#>  $ hat   :List of 1
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:99] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. .. ..@ p       : int [1:100] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. .. ..@ Dim     : int [1:2] 99 99
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:99] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..- attr(*, "dim")= int [1:2] 1 1
#>  $ params: num [1, 1:9] 99 0.2 1 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:9] "n" "mu" "sigma" "nonc" ...
```

We need a list of sparse matrices which
[`make_hat_mats()`](https://profoc.berrisch.biz/reference/make_hat_mats.md)
returns. So we can pass that directly to our class.

``` r
model$hat_pr <- hat$hat
```

The other smoothing matrices have to be filled with defaults (lists of
sparse identity matrices). Usually
[`online()`](https://profoc.berrisch.biz/reference/online.md) takes care
of this. But we can do it manually as well.

``` r
model$basis_mv <- list(Matrix::sparseMatrix(i = 1:D, j = 1:D, x = 1))
model$basis_pr <- list(Matrix::sparseMatrix(i = 1:P, j = 1:P, x = 1))
model$hat_mv <- list(Matrix::sparseMatrix(i = 1:D, j = 1:D, x = 1))
```

Now we can specify the parameter grid. We will stick to the defaults
here:

``` r
parametergrid <- as.matrix(
  expand.grid(
    forget_regret = 0,
    soft_threshold = -Inf,
    hard_threshold = -Inf,
    fixed_share = 0,
    basis_pr_idx = 1,
    basis_mv_idx = 1,
    hat_pr_idx = 1,
    hat_mv_idx = 1,
    gamma = 1,
    loss_share = 0,
    regret_share = 0
  )
)

model$params <- parametergrid
```

Finally, we can run `model$set_defaults()`. This populates initial
states (w0 for weights and R0 for regret).

``` r
model$set_defaults()
```

Now `model$set_grid_objects()` will create the grid objects
(performance, weights, regret etc.)

``` r
model$set_grid_objects()
```

Finally, we can run `model$learn()` to start the learning process.

``` r
model$learn()
```

## Accessing the results

The learning process fills the class objects. So we can inspect them
using the `$` operator, like we would with any other R object. For
example, we can access the weights:

``` r
head(model$weights[[T]][, , 1])
#> [1] 0.4892732 0.4058720 0.5979021 0.5770314 0.5487374 0.6287706
```

However, we can also use the post processing function of
[`online()`](https://profoc.berrisch.biz/reference/online.md) to access
the results. This will create output that is identical to the output of
[`online()`](https://profoc.berrisch.biz/reference/online.md):

``` r
names <- list(y = dimnames(y))
names$experts <- list(
  1:T,
  paste("Marginal", 1:D),
  tau,
  paste("Expert", 1:N)
)

output <- post_process_model(model, names)
```

We can now use `output` in
[`update()`](https://rdrr.io/r/stats/update.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and others.

At this point, we do not need to keep the model in memory anymore. So we
can delete it:

``` r
rm(model)
```

## Summary

The C++ class `conline` allows you to gain fine grained control over the
learning process. However, it is not as convenient as the
[`online()`](https://profoc.berrisch.biz/reference/online.md) wrapper.
So you should only use it if you need to alter the default behavior.
However, mixing up helper functions from
[`online()`](https://profoc.berrisch.biz/reference/online.md) and the
C++ class is possible. So you can compute your combinations using the
class interface while still being able to use
[`update()`](https://rdrr.io/r/stats/update.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and others
afterward.
