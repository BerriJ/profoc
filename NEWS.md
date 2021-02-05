profoc v0.5.X (Release date: 2021-XX-XX)
==============

Changes:

* The pinball_loss and loss_pred functions were replaced by a more flexible function called 'loss'.
* The weights object is changed from a ( T+1, K, P ) array to a ( T+1, P, K ) array to match other objects' dimensions. Now the following indexing scheme is consistent throughout the package: (Time, Probabilities, Experts, Parameter combination)
* Fixed Bug that caused single quantile calculations to fail.

profoc v0.5.2 (Release date: 2021-01-29)
==============

Changes:

* Initial release on GitHub