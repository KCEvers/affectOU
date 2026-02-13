# Plot OU simulation

Visualise an Ornstein-Uhlenbeck affect simulation using different types
of plots. See specific plotting functions for allowed arguments and
details.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
plot(x, type = c("time", "histogram", "acf", "phase"), ...)
```

## Arguments

- x:

  A `simulate_affectOU` model object produced by
  [`simulate.affectOU()`](https://kcevers.github.io/affectOU/reference/simulate.affectOU.md)

- type:

  Type of plot; one of `"time"`, `"histogram"`, `"acf"`, or `"phase"`

- ...:

  Additional parameters passed to specific plotting functions

## Value

NULL (invisibly), called for side effects only

## Details

Available plot types:

- `"time"`: Time series trajectories of affect dimensions produced by
  [`ou_plot_time()`](https://kcevers.github.io/affectOU/reference/ou_plot_time.md).

- `"histogram"`: Histograms of affect distributions produced by
  [`ou_plot_histogram()`](https://kcevers.github.io/affectOU/reference/ou_plot_histogram.md).

- `"acf"`: Autocorrelation and cross-correlation functions produced by
  [`ou_plot_acf()`](https://kcevers.github.io/affectOU/reference/ou_plot_acf.md).

- `"phase"`: Phase portraits produced by
  [`ou_plot_phase()`](https://kcevers.github.io/affectOU/reference/ou_plot_phase.md).

## Examples

``` r
# Simulate a 2-dimensional OU affect process
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 3)

# Plot simulation
plot(sim, type = "time")

plot(sim, type = "histogram")

plot(sim, type = "acf")

plot(sim, type = "phase")

```
