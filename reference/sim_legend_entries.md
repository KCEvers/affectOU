# Build legend entries for multiple simulations

Uses the lty assignment from
[`assign_sim_lty()`](https://kcevers.github.io/affectOU/reference/assign_sim_lty.md)
to produce one legend entry per occupied line-type group, labelled with
the simulation range (e.g. "Sim 1-4", "Sim 5-8").

## Usage

``` r
sim_legend_entries(nsim, lty_vec, col_sim, lwd = 2, sim_ids = seq_len(nsim))
```

## Arguments

- nsim:

  Number of simulations.

- lty_vec:

  Integer vector of length nsim from
  [`assign_sim_lty()`](https://kcevers.github.io/affectOU/reference/assign_sim_lty.md).

- col_sim:

  Character vector of colours (length nsim) from
  [`generate_shades()`](https://kcevers.github.io/affectOU/reference/generate_shades.md).

- lwd:

  Line width for legend entries.

## Value

Named list with `text`, `lty`, `lwd`, `col`, or `NULL` when `nsim <= 1`.
