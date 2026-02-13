# Print summary of affectOU model

Print summary of affectOU model

## Usage

``` r
# S3 method for class 'summary_affectOU'
print(x, digits = 3, max_dim = 20, ...)
```

## Arguments

- x:

  An object of class `summary_affectOU`

- digits:

  Number of digits to display

- max_dim:

  Maximum number of dimensions to display details for

- ...:

  Additional arguments (unused)

## Examples

``` r
model <- affectOU(ndim = 2)
summary_model <- summary(model)
print(summary_model)
#> 
#> ── 2D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: [0, 0]
#> SD: [1, 1]
#> Half-life: [1.386, 1.386]
#> Relaxation time (τ): [2, 2]
#> 
#> ── Structure ──
#> 
#> Coupling: none
#> Noise: independent
```
