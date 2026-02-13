# Plot observed vs fitted trajectory

Compare the observed data with the fitted trajectory over time. This
plot helps assess how well the model captures the overall pattern of the
data.

## Usage

``` r
ou_plot_fit_time(
  fit,
  palette = "Temps",
  alpha = 0.5,
  main = "Observed and Fitted Trajectory",
  xlab = "Time",
  ylab = "Affect",
  legend_position = "topright",
  ...
)
```

## Arguments

- fit:

  A `fit_affectOU` object

- palette:

  Color palette. Should be one
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html).

- alpha:

  Alpha transparency for colors (0 = transparent, 1 = opaque)

- main:

  Main title

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

NULL (invisibly), called for side effects only.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
ou_plot_fit_time(fitted)
```
