# Confidence intervals for fitted OU affect model

Confidence intervals for fitted OU affect model

## Usage

``` r
# S3 method for class 'fit_affectOU'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class
  [`fit_affectOU`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md).

- parm:

  Optional character vector of parameter names to include. If missing,
  all parameters are included.

- level:

  Confidence level for intervals (default 0.95).

- ...:

  Additional arguments (unused).

## Value

Matrix of confidence intervals with columns for lower and upper bounds.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
confint(fitted)
#>              2.5%     97.5%
#> theta  0.39396370 0.5704666
#> mu    -0.08125462 0.2786512
#> gamma  0.97002652 1.0097711
```
