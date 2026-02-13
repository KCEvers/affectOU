# Print relaxation time

Print relaxation time

## Usage

``` r
# S3 method for class 'relaxation_affectOU'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `relaxation_affectOU` object from
  [relaxation()](https://kcevers.github.io/affectOU/reference/relaxation.affectOU.md)

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

Returns `x` invisibly.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
print(relaxation(model))
#> 
#> ── Relaxation time of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Half-life (t₁/₂): 1.386 time units
#> Relaxation time (τ): 2 time units
#> 
#> Time for perturbation to decay to:
#>      50%     37%     14%      5%      1%
#>    1.386   2.000   4.000   6.000  10.000
```
