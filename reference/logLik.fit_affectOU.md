# Extract log-likelihood from fitted OU affect model

Extract log-likelihood from fitted OU affect model

## Usage

``` r
# S3 method for class 'fit_affectOU'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `fit_affectOU`

- ...:

  Additional arguments (unused)

## Value

Numeric log-likelihood value

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
logLik(fitted)
#> [1] -1206.394
```
