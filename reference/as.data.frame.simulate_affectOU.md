# Convert simulation results to data frame

Converts an `simulate_affectOU` object to a data frame in either long or
wide format.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  direction = c("long", "wide"),
  ...
)
```

## Arguments

- x:

  A `simulate_affectOU` object

- row.names:

  NULL or a character vector giving the row names for the data frame.
  Missing values are not allowed.

- optional:

  Logical. If TRUE, setting row names and converting column names is
  optional. Included for compatibility with the generic.

- direction:

  Character, either `"long"` (default) or `"wide"`. Long format has one
  row per time-dimension-simulation combination. Wide format has one row
  per time point with dimensions and simulations spread across columns.

- ...:

  Additional arguments (unused)

## Value

For `direction = "long"`, a data frame with columns:

- time:

  Observation time

- dim:

  Dimension index

- sim:

  Simulation index

- value:

  Simulated value

For `direction = "wide"`, a data frame with `time` as the first column
and subsequent columns named `dim{d}.sim{s}` (or simplified names for
univariate or single simulations).

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Long format (default) - one row per timepoint, dimension, and simulation
df_long <- as.data.frame(sim)
head(df_long)
#>   time dim sim       value
#> 1 0.00   1   1  0.00000000
#> 2 0.01   1   1 -0.12632048
#> 3 0.02   1   1 -0.24281192
#> 4 0.03   1   1 -0.23046567
#> 5 0.04   1   1 -0.04570529
#> 6 0.05   1   1 -0.15872705

# Wide format - one row per time point
df_wide <- as.data.frame(sim, direction = "wide")
head(df_wide)
#>   time   dim1.sim1   dim2.sim1   dim1.sim2 dim2.sim2   dim1.sim3  dim2.sim3
#> 1 0.00  0.00000000 0.000000000  0.00000000 0.0000000  0.00000000 0.00000000
#> 2 0.01 -0.12632048 0.001993167 -0.12195202 0.1362198 -0.03307022 0.16951448
#> 3 0.02 -0.24281192 0.048566532 -0.13139151 0.2058622 -0.03429512 0.09540621
#> 4 0.03 -0.23046567 0.033155678  0.01353667 0.2416916 -0.19353501 0.15481470
#> 5 0.04 -0.04570529 0.201903558  0.19025801 0.3814233 -0.11913414 0.11965166
#> 6 0.05 -0.15872705 0.052531116  0.27615935 0.2348255 -0.17774534 0.24242733
```
