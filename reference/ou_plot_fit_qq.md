# QQ plot of residuals

Visualize a normal Q-Q plot of the residuals to assess their
distributional properties. Ideally, residuals should lie approximately
along the reference line, indicating normality.

## Usage

``` r
ou_plot_fit_qq(
  fit,
  palette = "Temps",
  main = "Normal Q-Q Plot of Residuals",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles",
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
ou_plot_fit_qq(fitted)
```
