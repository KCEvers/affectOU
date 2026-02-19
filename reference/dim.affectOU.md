# Extract dimensionality

Extract the number of dimensions of an affect Ornstein-Uhlenbeck (OU)
model.

## Usage

``` r
# S3 method for class 'affectOU'
dim(x)
```

## Arguments

- x:

  An object of class `affectOU`

## Value

Integer, the dimensionality of the process.

## Examples

``` r
model <- affectOU()
dim(model)
#> [1] 1
```
