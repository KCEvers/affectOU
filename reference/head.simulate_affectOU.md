# Head of simulation results

Returns the first `n` time points from a `simulate_affectOU` object.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
head(x, n = 6L, ...)
```

## Arguments

- x:

  A `simulate_affectOU` simulation object

- n:

  Number of time points to keep from the start (default 6)

- ...:

  Additional arguments passed to
  [`utils::head()`](https://rdrr.io/r/utils/head.html)

## Value

A simulation data.frame truncated to the first `n` time points.

## Examples

``` r
model <- affectOU()
sim <- simulate(model)
head(sim)
#>   time dim sim      value
#> 1 0.00   1   1 -0.7684690
#> 2 0.01   1   1 -0.8854581
#> 3 0.02   1   1 -0.8332532
#> 4 0.03   1   1 -0.8853885
#> 5 0.04   1   1 -0.9527939
#> 6 0.05   1   1 -0.9726395
```
