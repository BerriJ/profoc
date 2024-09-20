profoc 1.3.3
==============

## Improvements

* We adjusted the integration of `rcpptimer`. This simplifies the code and makes use of the API of `rcpptimer` 1.2.0 which is expected to be stable. 

profoc 1.3.2
==============

## Improvements

* The timer functionality of online was moved to a seperate package [rcpptimer](https://github.com/BerriJ/rcpptimer). This is now added to profoc as a dependency. The timer-related code was removed. This makes the code more modular and easier to maintain. The timer functionality is now also available for other 'R' packages and even other languages (i.e. Python, via [cpptimer](https://github.com/BerriJ/cpptimer) and [cppytimer](https://github.com/BerriJ/cppytimer)).

profoc 1.3.1
==============

## Improvements
* Adjusted the clock.h code so that a larger share of code can be shared between the R and Python versions of that file.
* clock.h now uses welfords online algorithm to calculate the mean and variance of the timings. SD is reported in the times table.

## Fixes
* Fixed an integer overflow in the clock.h code which caused the package to fail on some systems.
* Fixed online() function for cases where the regret is exactly zero. This can happen if:
* * Only a single expert is used
* * Only two experts are provided and they both have the same predictions (in the beginning).

profoc 1.3.0
==============

## Improvements

* New articles explain how to use methods on `online()` objects to deploy online learning algorithms in production.
* The `conline` C++ class now exposes `weights` to R.
* A new article on the usage of the `conline` C++ class was added.
* Various functions are now exported to R to allow easier interaction with the `conline` C++ class. These functions are: `init_experts_list()`, `make_basis_mats` and `make_hat_mats`
* The code of `online()` was simplified a bit by utilizing the new `init_experts_list()` function.
* Function `post_process_model()` was improved and is now exposed to be used in conjunction with the `conline` C++ class.
* Move aggregation of timings from cppclock.R to clock.h. This make it faster, easier to maintain and simplifies the code (which will be used in python in the future as well).

profoc 1.2.1
==============

## Improvements

* `online()` outputs now include `predictions_got_sorted`. A matrix which indicates whether quantile crossing occured and predictions have been sorted.
* `tidy()` methods were added to convert `weights`, `predictions` and loss objects of `online()` output to a tibble (for further analysis, plotting etc.)
* A [Get started](https://profoc.berrisch.biz/articles/profoc.html) article was added to the docs.
* Docs of the development version were added to the [website](https://profoc.berrisch.biz/dev/)

## Fixes

* This release fixes import / export of of the `autoplot()` method. In consequence, ggplot2 became a new dependency of this package.

profoc 1.2.0
==============

## Improvements:

* Periodic splines and penalties added for smoothing the weights in `online()`.

### Internal changes
* `profoc` now depends on `R >= 4.3.0` to ensure C++17 support.

profoc 1.1.1
==============

## Fixes:

* Distribution of the knots is now correct for `ncp < 0`. 

profoc 1.1.0
==============

## Improvements:

* New `penalty()` function which works with equidistant and non-equidistant knots. 

## Fixes:

* Calculation of the P-Spline penalty if non-uniform B-Splines are used.

profoc 1.0.0
==============

## Changes:
* Now, `online()` saves memory by not reporting `past_performance` and `past_predictions_grid`. However, the cumulative performance and the most recent predictions w.r.t to the parameter grid are always included in the output. The former is used internally for choosing the best hyperparameter set, and the latter for updating the weights. Depending on the data and the parameter space considered, both objects may get large. You can still opt-in to include them in the output by setting `save_past_performance = TRUE` and `save_past_predictions_grid = TRUE` in `online()`.

### Internal changes
* Minor fixes and improvements to `online()` to reduce memory usage.

profoc 0.9.5
==============

## Internal changes
* Now `online()` is able to sample from grids of up to 2^64-1 rows.
* The new cpp sampling function `sample_int()` works similar to `sample.int()` and also respects seeds set by `set.seed()`. 

profoc 0.9.4
==============

## Fixes:

* `parametergrids` lets you provide custom grids of parameters in `online()` 

### Internal changes
* Significantly improved the initialization efficiency in `online()` when using large grids of parameters

profoc 0.9.3
==============

## Fixes:

* `forget_past_performance` had no effect in `online()` 
* Improved and fixed documentation

profoc 0.9.2
==============

## Fixes:

* Resolved a problem with plotting multivariate probabilistic models

profoc 0.9.1
==============

## Changes:

* Basis matrices are created differently now. This solves an issue where basis functions did not always sum to 1 when non-equidistant knot sequences were used.

profoc 0.9.0
==============

## Changes:

* `online` can now be used with multivariate data
  * Just pass a TxK matrix as `y` and a TxDxPxK array as `experts`
* Smoothing was improved. See the documentation for details on the revised interface.
* `summary.online` can be used to obtain selected parameters of `online` models

### Internal changes

* `online` uses Rcpp Modules to bundle data and functionality into an exposed C++ class
* Improvements to plot methods

profoc 0.8.5
==============

## Changes:

* `initial_weights` argument is replaced by `init`
  * `init` takes a named list and currently `intial_weights` and `R0` the initial weights and the initial cumulative regret can be provided. They have to be PxK or 1xK. 

### Internal changes

* Internal changes to improve readability
* Resolve C++ compilation warnings

profoc 0.8.4
==============

* Remove unused C++17 dependency

profoc 0.8.3
==============

## Changes:

* The `profoc` function was extended:
  * `regret` can now be passed as an array as before, or as a list, e.g. `list(regret = regret_array, share = 0.2)` if the provided regret should be mixed with the regret calculated by online.
  * `loss` can also be provided as a list, see above.

* The `batch` function can now minimize an alternative objection function using the quantile weighted CRPS
  * This weighting scheme can be activated by setting `qw_crps=TRUE`
  * It defaults to FALSE due to better performance
  
### Internal changes

* Changes to the splines2 package required some changes on our side. This affected the creation of the b-spline basis. See https://github.com/BerriJ/profoc/pull/3

profoc 0.8.0
==============

## Changes:

* First release on CRAN
* The `profoc` function was renamed to `online` for consistency.
* Added `batch` function to apply batch-learning.
* Added `oracle` function to approximate the oracle.
* Update, predict, plot, and print methods were added for `online` and `batch` objects.

### Interface:

* Unfortunately, we decided to apply significant changes to the API. This likely breaks old code. Please refer to the respective function documentation for more information.

### Internal changes:

* The b-spline basis is now calculated using a fast C++ function imported from the [splines2](https://github.com/wenjie2wang/splines2) R package.
* The source code is now distributed across different files.

profoc 0.7.0
==============

## Changes:

The spline functions where rewritten to add the ability of using a non-equidistant knot sequence and a penalty term defined on the Sobolev space. This change induces breaking changes to small parts of the API.

### Interface:

* `ndiff` defines the degree of differencing for creating the penalty term. For values between 1 and 2 a weighted sum of the difference penalization matrices is used.
* `rel_nseg` is replaced by `knot_distance` ( distance between knots). Defaults to 0.025, which corresponds to the grid steps when knot_distance_power = 1 (the default).
* A new parameter `knot_distance_power` defines if knots are uniformly distributed. Defaults to 1, which corresponds to the equidistant case. Values less than 1 create more knots in the center, while values above 1 concentrate more knots in the tails.
* A new parameter `allow_quantile_crossing` defines if quantile crossing is allowed. Defaults to false, which means that predictions will be sorted.

### Internal changes:

* Functions for calculating the b-spline basis are now exported to R as internal functions of the package. They can be accessed using the `package:::function` notation.

profoc 0.6.0
==============

## Changes:

### Interface:

* `y` must now be a matrix of either $\text{T} \times 1$ or $\text{T} \times \text{P}$. 
* `trace` specifies whether a progress bar will be printed or not. Default to `TRUE`.
* `loss_function` lets you now specify "quantile", "expectile" or "percentage". All functions are generalized as in [Gneitling 2009](https://arxiv.org/abs/0912.0902). The power can be scaled by `loss_parameter`. The latter defaults to 1, which leads to the well-known quantile, squared, and absolute percentage loss.
* `gradient` lets you specify whether the learning algorithm should consider actual loss or a linearized version using the gradient of the loss. Defaults to `TRUE` (gradient-based learning).
* `forget_performance` was added. It defines the share of the past performance that will be ignored when selecting the best parameter combination.
* Renamed `forget` parameter to `forget_regret` to underline its reference to the regret.
* New `init_weights` parameter. It has to be either a Kx1 or KxP matrix specifying the experts' starting weights.
* Add `lead_time` parameter. offset for expert forecasts. Defaults to 0 which means that experts predict t+1 at t. Setting this to h means experts' predictions refer to t+1+h at time t. The weight updates delay accordingly.

### Internal changes:

* If more expert forecasts than observations are provided, the excessive expert forecasts are used for prediction using the most recent weights.
* `tau` is now optional. It defaults to 1:P/(P+1). A scalar given to tau will be repeated P times. The latter is useful in multivariate settings.
* The `pinball_loss` and `loss_pred` functions were replaced by a more flexible function called `loss`.
* The `weights` object is changed from a $(\text{T}+1 \times \text{K} \times \text{P})$ array to a $(\text{T}+1 \times \text{P} \times \text{K})$ array to match other objects' dimensions. Now the following indexing scheme is consistent throughout the package: (Time, Probabilities, Experts, Parameter combination)
* Fixed Bug that caused single quantile calculations to fail.
* Various internal changes to improve readability and performance.

profoc 0.5.2
==============

## Changes:

* Initial release on GitHub