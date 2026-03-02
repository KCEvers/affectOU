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
#>   time dim sim      value
#> 1 0.00   1   1 -0.5914377
#> 2 0.01   1   1 -0.6800363
#> 3 0.02   1   1 -0.7285889
#> 4 0.03   1   1 -0.7447949
#> 5 0.04   1   1 -0.7889493
#> 6 0.05   1   1 -0.7532798

# Wide format - one row per time point
df_wide <- as.data.frame(sim, direction = "wide")
head(df_wide)
#>   time  dim1.sim1 dim2.sim1  dim1.sim2  dim2.sim2     dim1.sim3  dim2.sim3
#> 1 0.00 -0.5914377 0.5524307 -0.9636726 -0.7148733  0.0458320665 -0.6324440
#> 2 0.01 -0.6800363 0.5015683 -1.0327679 -0.9399898 -0.0189821663 -0.4979532
#> 3 0.02 -0.7285889 0.4336205 -0.8400477 -0.9135680  0.0869041476 -0.5740382
#> 4 0.03 -0.7447949 0.4865931 -0.7611739 -1.0446130 -0.0120263762 -0.7609862
#> 5 0.04 -0.7889493 0.5326139 -0.8675570 -0.8617633 -0.2065835213 -0.7989514
#> 6 0.05 -0.7532798 0.3339744 -0.8590048 -0.9158483  0.0004971437 -0.7561335
```
