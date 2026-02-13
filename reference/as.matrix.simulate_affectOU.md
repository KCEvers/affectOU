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
#>      time dim sim       value
#> [1,] 0.00   1   1  0.00000000
#> [2,] 0.01   1   1  0.02924642
#> [3,] 0.02   1   1  0.03342802
#> [4,] 0.03   1   1  0.04899601
#> [5,] 0.04   1   1 -0.15165722
#> [6,] 0.05   1   1 -0.17943435

# Wide format
mat_wide <- as.matrix(sim, direction = "wide")
head(mat_wide)
#>        dim1.sim1  dim2.sim1   dim1.sim2   dim2.sim2  dim1.sim3    dim2.sim3
#> 0     0.00000000 0.00000000  0.00000000  0.00000000  0.0000000 0.0000000000
#> 0.01  0.02924642 0.01470536 -0.08944247 -0.07456070 -0.1104323 0.0002603787
#> 0.02  0.03342802 0.17518862  0.12337504 -0.14724970 -0.1899574 0.0335518376
#> 0.03  0.04899601 0.42753036  0.09238958  0.05033636 -0.3087358 0.1488129498
#> 0.04 -0.15165722 0.54516828  0.12244803  0.17586282 -0.4818593 0.3050938455
#> 0.05 -0.17943435 0.55616230  0.08447718  0.24068578 -0.3883452 0.0599221413
```
