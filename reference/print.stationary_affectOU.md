# Print stationary distribution

Print stationary distribution

## Usage

``` r
# S3 method for class 'stationary_affectOU'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `stationary_affectOU` object from
  [stationary()](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md)

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

Returns `x` invisibly.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
print(stationary(model))
#> 
#> ── Stationary distribution of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Mean: 0
#> SD: 1
#> 95% interval: [-2, 2]
```
