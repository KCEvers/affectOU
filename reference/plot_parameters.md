# Plot parameter documentation

Plot parameter documentation

## Arguments

- which_dim:

  Dimension indices to plot (NULL for all)

- which_sim:

  Simulation indices to plot (NULL for all)

- by_dim:

  Logical; plot each dimension in separate panel?

- palette:

  Color palette. Should be one
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html).

- alpha:

  Alpha transparency for colors (0 = transparent, 1 = opaque)

- share_xaxis:

  Logical; use same x-axis limits for all panels?

- share_yaxis:

  Logical; use same y-axis limits for all panels?

- freq:

  Logical; plot frequency instead of density?

- breaks:

  Number of histogram breaks

- lag.max:

  Maximum lag to compute. Specified in terms of saved time points. For
  example, `lag.max = 10` corresponds to 10 time units and 100 lags with
  `save_at = 0.1`.

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
