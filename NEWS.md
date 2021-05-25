profoc 0.7.0 (Release date: 2021-05-XX)
==============

## Changes:

The spline functions where rewritten to add the ability of using a non-equidistant knot sequence and a penalty term defined on the sobolev space. This change induces breaking changes to small parts of the API.

### Interface:

* `ndiff` defines the degree of differencing for creating the penalty term. For values between 1 and 2 a weighted sum of the difference penalization matrices is used.
* `rel_nseg` is replaced by `knot_distance` ( distance between knots). Defaults to 0.025, which corresponds to the grid steps when knot_distance_power = 1 (the default).
* A new parameter `knot_distance_power` defines if knots are uniformly distributed. Defaults to 1 which corresponds to the equidistant case. Values less than 1 create more knots in the center while values above 1 concentrate more knots in the tails.

### Internal changes:

* Functions for calculating the b-spline basis are now exported to R as internal functions of the package. They can be accessed using the `package:::function` notation.

profoc 0.6.0 (Release date: 2021-03-04)
==============

## Changes:

### Interface:

* `y` must now be a matrix of either $\text{T} \times 1$ or $\text{T} \times \text{P}$. 
* `trace` specifies whether a progress bar will be printed or not. Default to `TRUE`.
* `loss_function` lets you now specify "quantile", "expectile" or "percentage". All functions are generalized as in [Gneitling 2009](https://arxiv.org/abs/0912.0902). The power can be scaled by `loss_parameter`. The latter defaults to 1, which leads to the well-known quantile, squared, and absolute percentage loss.
* `gradient` lets you specify wether the learning alrogithms should consider actual loss or a linearized version using the gradient of the loss. Defaults to `TRUE` (gradient based learning).
* `forget_performance` was added. It defines the share of the past performance that will be ignored when selecting the best parameter combination.
* Renamed `forget` parameter to `forget_regret` to underline its reference to the regret.
* New `init_weights` parameter. Has to be either a Kx1 or KxP matrix specifying the experts' starting weights.
* Add `lead_time` parameter. offset for expert forecasts. Defaults to 0 which means that experts predict t+1 at t. Setting this to h means experts predictions refer to t+1+h at time t. The weight updates delay accordingly.

### Internal changes:

* If more expert forecasts than observations are provided, the excessive expert forecasts are used for prediction using the most recent weights.
* `tau` is now optional. It defaults to 1:P/(P+1). A scalar given to tau will be repeated P times. The latter is useful in multivariate settings.
* The pinball_loss and loss_pred functions were replaced by a more flexible function called `loss`.
* The `weights` object is changed from a $(\text{T}+1 \times \text{K} \times \text{P})$ array to a $(\text{T}+1 \times \text{P} \times \text{K})$ array to match other objects' dimensions. Now the following indexing scheme is consistent throughout the package: (Time, Probabilities, Experts, Parameter combination)
* Fixed Bug that caused single quantile calculations to fail.
* Various internal changes to improve readability and performance.

profoc 0.5.2 (Release date: 2021-01-29)
==============

## Changes:

* Initial release on GitHub