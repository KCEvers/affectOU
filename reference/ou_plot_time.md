# Plot simulation trajectory

Visualise the time series trajectories of affect dimensions from an OU
affect simulation. Each dimension is plotted in a separate panel, with
multiple simulations overlaid within each panel.

## Usage

``` r
ou_plot_time(
  x,
  which_dim = NULL,
  which_sim = NULL,
  by_dim = TRUE,
  palette = "Temps",
  col_theory = "grey30",
  alpha = 1,
  share_yaxis = TRUE,
  main = "Affect Dynamics",
  sub = paste("Dimension", if (is.null(which_dim)) {
    
    seq.int(x[["model"]][["ndim"]])
 } else {
     which_dim
 }),
  xlab = "Time",
  ylab = "Affect",
  legend_position = "topright",
  ...
)
```

## Arguments

- x:

  A `simulate_affectOU` model object produced by
  [`simulate.affectOU()`](https://kcevers.github.io/affectOU/reference/simulate.affectOU.md)

- which_dim:

  Dimension indices to plot (NULL for all)

- which_sim:

  Simulation indices to plot (NULL for all)

- by_dim:

  Logical; plot each dimension in separate panel?

- palette:

  Color palette. Should be one
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html).

- col_theory:

  Color for `mu` (i.e., attractor) line

- alpha:

  Alpha transparency for colors (0 = transparent, 1 = opaque)

- share_yaxis:

  Logical; use same y-axis limits for all panels?

- main:

  Main title

- sub:

  Subtitle for panels

- xlab:

  X-axis label

- ylab:

  Y-axis label

- legend_position:

  Position of legend (one of `"bottomright"`, `"bottom"`,
  `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`,
  `"right"`, `"center"`)

- ...:

  Additional graphical parameters

## Value

NULL (invisibly), called for side effects only

## Attractor Line

The horizontal dashed line shows the attractor level \\\mu\\.
Trajectories fluctuate around this baseline, pulled back by the drift
term \\\theta(\mu - X(t))\\. The strength of mean reversion (\\\theta\\)
determines how tightly trajectories cluster around \\\mu\\.

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)
ou_plot_time(sim)


# Plot dimensions in one panel
sim <- simulate(model, nsim = 1)
ou_plot_time(sim, by_dim = FALSE)

```
