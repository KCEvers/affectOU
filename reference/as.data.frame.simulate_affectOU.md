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
and subsequent columns named `dim{d}` when `nsim = 1`, or
`dim{d}.sim{s}` when `nsim > 1`.

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Long format (default) - one row per timepoint, dimension, and simulation
df_long <- as.data.frame(sim)
head(df_long)
#>   time dim sim     value
#> 1 0.00   1   1 -1.132503
#> 2 0.01   1   1 -1.119858
#> 3 0.02   1   1 -1.059016
#> 4 0.03   1   1 -1.101821
#> 5 0.04   1   1 -1.161752
#> 6 0.05   1   1 -1.100803

# Wide format - one row per time point
df_wide <- as.data.frame(sim, direction = "wide")
head(df_wide)
#>   time dim1.sim1 dim2.sim1  dim1.sim2   dim2.sim2  dim1.sim3  dim2.sim3
#> 1 0.00 -1.132503 -1.483629 -0.9129529 -0.15528623 -0.4677025 -0.8242106
#> 2 0.01 -1.119858 -1.535355 -0.9212173 -0.25087706 -0.4063142 -0.8155064
#> 3 0.02 -1.059016 -1.619234 -0.9880985 -0.32353629 -0.4675271 -0.8760139
#> 4 0.03 -1.101821 -1.663090 -1.2118489 -0.13436227 -0.3338609 -0.7658424
#> 5 0.04 -1.161752 -1.674624 -1.1840678 -0.05901687 -0.4107663 -0.8605092
#> 6 0.05 -1.100803 -1.714129 -1.3137603 -0.16891078 -0.5985307 -1.0508240
```
