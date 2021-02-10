profoc v0.5.X (Release date: 2021-XX-XX)
==============

Changes:

* The `init_weights` parameter has to be either a Kx1 or Kx99 matrix specifying the experts' starting weights.
* A new `trace` parameter lets you decide whether a progress bar is printed.
* `tau` is now optional. It defaults to 1:P/(P+1). A scalar given to tau will be repeated P times. The latter is useful in multivariate settings.
* `y` must now be a matrix of either $\text{T} \times 1$ or $\text{T} \times \text{P}$. 
* The pinball_loss and loss_pred functions were replaced by a more flexible function called `loss`.
* The `weights` object is changed from a $(\text{T}+1 \times \text{K} \times \text{P})$ array to a $(\text{T}+1 \times \text{P} \times \text{K})$ array to match other objects' dimensions. Now the following indexing scheme is consistent throughout the package: (Time, Probabilities, Experts, Parameter combination)
* Fixed Bug that caused single quantile calculations to fail.
* Various internal changes to improve the readability and performance of the code

profoc v0.5.2 (Release date: 2021-01-29)
==============

Changes:

* Initial release on GitHub