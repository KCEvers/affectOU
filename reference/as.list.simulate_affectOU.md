# Convert simulation results to list

Returns the simulation data as a list with one data frame per
simulation.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
as.list(x, direction = c("wide", "long"), ...)
```

## Arguments

- x:

  A `simulate_affectOU` object

- direction:

  Character string specifying the output format: `"wide"` (default)
  returns one column per dimension, `"long"` returns a single value
  column with a dimension indicator.

- ...:

  Additional arguments (unused)

## Value

A list of data frames, one per simulation. Each data frame contains:

- time:

  Time points

- dim1, dim2, ...:

  (wide format) Values for each dimension

- dim:

  (long format) Dimension indicator

- value:

  (long format) Simulated values

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Wide format (default): one column per dimension
lst_wide <- as.list(sim)
head(lst_wide[[1]])
#>   time       dim1      dim2
#> 1 0.00 -0.2036899 -1.182427
#> 2 0.01 -0.1850442 -1.127589
#> 3 0.02 -0.1028183 -1.230978
#> 4 0.03 -0.1129939 -1.159538
#> 5 0.04 -0.1981627 -1.017135
#> 6 0.05 -0.1529370 -1.072746

# Long format: single value column with dimension indicator
lst_long <- as.list(sim, direction = "long")
head(lst_long[[1]])
#>   time  dim      value
#> 1 0.00 dim1 -0.2036899
#> 2 0.01 dim1 -0.1850442
#> 3 0.02 dim1 -0.1028183
#> 4 0.03 dim1 -0.1129939
#> 5 0.04 dim1 -0.1981627
#> 6 0.05 dim1 -0.1529370
```
