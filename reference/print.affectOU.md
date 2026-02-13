# Print affect OU model

Print affect OU model

## Usage

``` r
# S3 method for class 'affectOU'
print(x, digits = 3, max_dim = 20, ...)
```

## Arguments

- x:

  An object of class `affectOU`

- digits:

  Number of digits to display

- max_dim:

  Maximum number of dimensions to display details for

- ...:

  Additional arguments (unused)

## Examples

``` r
model <- affectOU(ndim = 2)
print(model)
#> 
#> ── 2D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> dX(t) = Θ(μ − X(t))dt + Γ dW(t)
#> 
#> μ = [0.000, 0.000]
#> 
#> Θ:
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> Γ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> Σ = ΓΓᵀ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
