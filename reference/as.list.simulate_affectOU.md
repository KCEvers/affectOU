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
#>   time      dim1       dim2
#> 1 0.00 0.8130071 -1.0902676
#> 2 0.01 0.7982523 -1.0195320
#> 3 0.02 0.7085274 -0.8778291
#> 4 0.03 0.7492196 -0.9341361
#> 5 0.04 0.6924288 -0.9359415
#> 6 0.05 0.6721075 -0.9771002

# Long format: single value column with dimension indicator
lst_long <- as.list(sim, direction = "long")
head(lst_long[[1]])
#>   time  dim     value
#> 1 0.00 dim1 0.8130071
#> 2 0.01 dim1 0.7982523
#> 3 0.02 dim1 0.7085274
#> 4 0.03 dim1 0.7492196
#> 5 0.04 dim1 0.6924288
#> 6 0.05 dim1 0.6721075
```
