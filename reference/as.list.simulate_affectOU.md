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
#>   time         dim1         dim2
#> 1 0.00  0.000000000  0.000000000
#> 2 0.01  0.113245781  0.002042356
#> 3 0.02  0.128889179  0.040807355
#> 4 0.03 -0.006977228 -0.098792882
#> 5 0.04 -0.124309674 -0.146183066
#> 6 0.05 -0.136735460 -0.118359125

# Long format: single value column with dimension indicator
lst_long <- as.list(sim, direction = "long")
head(lst_long[[1]])
#>   time  dim        value
#> 1 0.00 dim1  0.000000000
#> 2 0.01 dim1  0.113245781
#> 3 0.02 dim1  0.128889179
#> 4 0.03 dim1 -0.006977228
#> 5 0.04 dim1 -0.124309674
#> 6 0.05 dim1 -0.136735460
```
