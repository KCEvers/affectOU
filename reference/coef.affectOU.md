# Extract model parameters

Extract the parameters of an affect Ornstein-Uhlenbeck (OU) model as a
list.

## Usage

``` r
# S3 method for class 'affectOU'
coef(object, ...)
```

## Arguments

- object:

  An object of class `affectOU`.

- ...:

  Additional arguments (unused).

## Value

A list containing the model parameters: `theta`, `mu`, `gamma`, and
`sigma`. For 1D models, these are returned as numeric scalars. For
multivariate models, they are returned as matrices.

## See also

[`affectOU()`](https://kcevers.github.io/affectOU/reference/affectOU.md)
to create an OU model and
[update()](https://kcevers.github.io/affectOU/reference/update.affectOU.md)
to modify model parameters.

## Examples

``` r
model <- affectOU(ndim = 2)
coef(model)
#> $theta
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> $mu
#> [1] 0 0
#> 
#> $gamma
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> $sigma
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
```
