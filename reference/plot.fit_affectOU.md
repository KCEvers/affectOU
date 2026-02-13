# Plot fitted model diagnostics

Provides various diagnostic plots for fitted
[`fit_affectOU`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md)
objects, including:

- `"time"`: Observed vs fitted trajectory over time produced by
  [`ou_plot_fit_time()`](https://kcevers.github.io/affectOU/reference/ou_plot_fit_time.md).

- `"residuals"`: Residuals (data - fitted) over time produced by
  [`ou_plot_fit_residuals()`](https://kcevers.github.io/affectOU/reference/ou_plot_fit_residuals.md).

- `"acf"`: ACF of residuals produced by
  [`ou_plot_fit_acf()`](https://kcevers.github.io/affectOU/reference/ou_plot_fit_acf.md).

- `"qq"`: Normal Q-Q plot of residuals produced by
  [`ou_plot_fit_qq()`](https://kcevers.github.io/affectOU/reference/ou_plot_fit_qq.md).

## Usage

``` r
# S3 method for class 'fit_affectOU'
plot(x, type = c("time", "residuals", "acf", "qq"), ...)
```

## Arguments

- x:

  An object of class
  [`fit_affectOU`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md)

- type:

  Type of plot; one of `"time"`, `"residuals"`, `"acf"`, or `"qq"`

- ...:

  Additional graphical parameters passed to the specific plotting
  functions

## Value

NULL (invisibly), called for side effects only.

## Examples

``` r
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.1)
data <- as.data.frame(sim)
fitted <- fit(model, data = data$value, times = data$time)
plot(fitted, type = "time")

plot(fitted, type = "residuals")

plot(fitted, type = "acf")

plot(fitted, type = "qq")
```
