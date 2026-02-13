# Extract model coefficients

Extract the parameters of an affect OU model as a list.

## Usage

``` r
# S3 method for class 'affectOU'
coef(object, ...)
```

## Arguments

- object:

  An object of class `affectOU`

- ...:

  Additional arguments (unused)

## Value

A list containing the model parameters: `theta`, `mu`, `gamma`, and
`sigma`. For 1D models, these are returned as numeric scalars. For
multivariate models, they are returned as matrices.

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
