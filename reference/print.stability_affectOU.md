# Print stability analysis

Print stability analysis

## Usage

``` r
# S3 method for class 'stability_affectOU'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `stability_affectOU` object from
  [stability()](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

Returns `x` invisibly.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
print(stability(model))
#> 
#> ── Stability analysis of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Stable (node). Deviations from the attractor decay exponentially.
```
