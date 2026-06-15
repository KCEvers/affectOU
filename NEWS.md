# affectOU 1.0.2

* Fixed `summary()` for multivariate simulations so means, standard deviations,
  covariances, and correlations are pooled by dimension correctly when multiple
  simulations are summarized.
* `simulate()` now always includes the requested `stop` time without returning
  `NA` values when `dt` and `save_at` do not divide the total simulation time
  exactly.
* `simulate()` now works for stable multivariate models with singular stationary
  covariance matrices, such as models with a dimension that has no stochastic
  noise.
* `affectOU()` now accepts singular positive semi-definite `sigma` matrices and
  gives a clear error for non-symmetric `sigma` matrices.
* `fit()` now gives a clear error when `times` are non-finite, duplicated, or
  not strictly increasing.
* `logLik()` for fitted models now returns a standard `logLik` object with
  degrees of freedom and number of observations attached.
* ACF plots for non-stationary simulations no longer show a theoretical
  stationary line.
* Corrected documentation for the multivariate stationary covariance equation
  and the `simpaper` example equation.

# affectOU 1.0.1

* `plot()` now supports `legend_position = "none"`

# affectOU 1.0.0

First stable version.
