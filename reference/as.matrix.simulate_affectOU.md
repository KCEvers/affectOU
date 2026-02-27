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

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Long format (default)
mat_long <- as.matrix(sim)
head(mat_long)
#>      time dim sim     value
#> [1,] 0.00   1   1 0.6756908
#> [2,] 0.01   1   1 0.6225849
#> [3,] 0.02   1   1 0.8554543
#> [4,] 0.03   1   1 0.7873109
#> [5,] 0.04   1   1 0.6592636
#> [6,] 0.05   1   1 0.7403810

# Wide format
mat_wide <- as.matrix(sim, direction = "wide")
head(mat_wide)
#>      dim1.sim1 dim2.sim1  dim1.sim2  dim2.sim2  dim1.sim3  dim2.sim3
#> 0    0.6756908 0.5082894 -0.6098916 -0.1268938 -0.6543554 -1.0839112
#> 0.01 0.6225849 0.6375182 -0.5868591 -0.1792499 -0.6350719 -1.1266001
#> 0.02 0.8554543 0.7115427 -0.3952107 -0.2649401 -0.5623341 -0.9640450
#> 0.03 0.7873109 0.7132281 -0.3928531 -0.1534357 -0.6168915 -0.9525534
#> 0.04 0.6592636 0.7132079 -0.3076045  0.1061817 -0.6697355 -1.0024722
#> 0.05 0.7403810 0.6309112 -0.3123168  0.1295897 -0.6829563 -0.9640234
```
