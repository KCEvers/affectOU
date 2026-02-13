# Calculate layout dimensions for multi-panel plots

Calculate layout dimensions for multi-panel plots

## Usage

``` r
get_layout(n, P, user_args = list(), by_dim = TRUE)
```

## Arguments

- n:

  Number of panels needed

- P:

  Parameter list containing nrow, ncol, mfrow, or mfcol

- user_args:

  Original user arguments (to check for mfrow/mfcol)

- by_dim:

  Logical; plot each dimension in separate panel?

## Value

List with nrow and ncol
