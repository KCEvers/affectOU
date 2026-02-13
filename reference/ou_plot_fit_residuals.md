# Plot residuals over time

Visualize the residuals (differences between observed and fitted values)
over time to assess model fit. Ideally, residuals should show no
systematic patterns and be randomly scattered around zero.

## Usage

``` r
ou_plot_fit_residuals(
  fit,
  palette = "Temps",
  main = "Residuals over Time",
  xlab = "Time",
  ylab = "Residual",
  alpha = 0.8,
  ...
)
```

## Arguments

- fit:

  A `fit_affectOU` object

- palette:

  Color palette. Should be one
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html).

- main:

  Main title

- xlab:

  X-axis label

- ylab:

  Y-axis label

- alpha:

  Alpha transparency for colors (0 = transparent, 1 = opaque)

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
ou_plot_fit_residuals(fitted)
```
