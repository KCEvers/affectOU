# Assign line types to simulations in contiguous groups

Divides `nsim` simulations into up to `n_lty` equal-sized contiguous
blocks and returns a vector of line-type indices (integers 1 to `n_lty`)
of length `nsim`. When `nsim <= n_lty` each simulation gets its own
unique line type.

## Usage

``` r
assign_sim_lty(nsim, n_lty = 5L)
```

## Arguments

- nsim:

  Number of simulations.

- n_lty:

  Number of distinct line types available (default 5).

## Value

Integer vector of length `nsim`.
