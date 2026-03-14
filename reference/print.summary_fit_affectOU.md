# Print summary of fitted OU affect model

Print summary of fitted OU affect model

## Usage

``` r
# S3 method for class 'summary_fit_affectOU'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class
  [`summary_fit_affectOU`](https://kcevers.github.io/affectOU/reference/summary.fit_affectOU.md).

- digits:

  Number of digits for numeric display.

- ...:

  Additional arguments (unused).

## Value

Returns `x` invisibly.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
print(summary(fitted))
#> 
#> ── Fitted Ornstein-Uhlenbeck Model Summary ─────────────────────────────────────
#> Method: mle, 5001 observations
#> 
#> ── Coefficients ──
#> 
#>       Estimate    SE          95% CI
#> theta    0.549 0.048  [0.455, 0.644]
#> mu      -0.011 0.082 [-0.172, 0.150]
#> gamma    1.010 0.010  [0.989, 1.030]
#> 
#> ── Goodness of fit ──
#> 
#> Log-likelihood: -1250.877
#> RMSE: 0.311
```
