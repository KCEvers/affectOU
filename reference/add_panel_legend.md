# Add legend inside a specific panel

Add legend inside a specific panel

## Usage

``` r
add_panel_legend(
  position = "topright",
  xlim,
  ylim,
  panel_row = 1,
  panel_col = 1,
  inset = c(0.02, 0.02),
  bg = "white",
  ...
)
```

## Arguments

- position:

  Legend position within the panel (e.g., "topright")

- xlim:

  X-axis limits for the target panel

- ylim:

  Y-axis limits for the target panel

- panel_row:

  Panel row index (1-based)

- panel_col:

  Panel column index (1-based)

- inset:

  Inset distance from margins as fraction of plot region

- bg:

  Background colour for the legend box

- ...:

  Additional arguments passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html)
