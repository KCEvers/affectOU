# Extract coefficients from fitted OU affect model

Extract coefficients from fitted OU affect model

## Usage

``` r
# S3 method for class 'fit_affectOU'
coef(object, ...)
```

## Arguments

- object:

  An object of class `fit_affectOU`

- ...:

  Additional arguments (unused)

## Value

Named vector of fitted parameters

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
coef(fitted)
#>      theta         mu      gamma 
#> 0.44684807 0.02906339 1.00871397 
```
