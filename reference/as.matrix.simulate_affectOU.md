# Convert simulation results to matrix

Returns the simulation data as a 2-dimensional matrix in either long or
wide format.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
as.matrix(x, direction = c("long", "wide"), ...)
```

## Arguments

- x:

  A `simulate_affectOU` object

- direction:

  Character, either `"long"` (default) or `"wide"`. Long format has one
  row per observation with columns for time, dim, sim, and value. Wide
  format has time as rows and dimension-simulation combinations as
  columns.

- ...:

  Additional arguments (unused)

## Value

For `direction = "long"`, a matrix with columns `time`, `dim`, `sim`,
`value`. For `direction = "wide"`, a matrix with dimensions (ntime x
ndim\*nsim). Columns iterate over dimensions first, then simulations.
Column names are `dim1`, `dim2`, ... when `nsim = 1`, or `dim1.sim1`,
`dim2.sim1`, ..., `dim1.sim2`, ... when `nsim > 1`.

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Long format (default)
mat_long <- as.matrix(sim)
head(mat_long)
#>      time dim sim       value
#> [1,] 0.00   1   1 -0.10486312
#> [2,] 0.01   1   1 -0.17606609
#> [3,] 0.02   1   1 -0.20910665
#> [4,] 0.03   1   1 -0.15723218
#> [5,] 0.04   1   1 -0.02467579
#> [6,] 0.05   1   1  0.05265964

# Wide format
mat_wide <- as.matrix(sim, direction = "wide")
head(mat_wide)
#>        dim1.sim1  dim2.sim1 dim1.sim2 dim2.sim2  dim1.sim3 dim2.sim3
#> 0    -0.10486312 -1.0508508 -1.350156 0.6629672 -0.2450390 0.5682976
#> 0.01 -0.17606609 -1.0276331 -1.456984 0.7318991 -0.3020935 0.4302807
#> 0.02 -0.20910665 -0.9549258 -1.416720 0.6672504 -0.1749767 0.3626938
#> 0.03 -0.15723218 -0.9998786 -1.422326 0.6838972 -0.2824929 0.3768920
#> 0.04 -0.02467579 -0.7588969 -1.468205 0.8691918 -0.3291890 0.4445699
#> 0.05  0.05265964 -0.8189686 -1.547450 0.8652275 -0.1706209 0.3849780
```
