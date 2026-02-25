
<!-- README.md is generated from README.Rmd. Please edit that file -->

# affectOU: Ornstein-Uhlenbeck model for Affect Dynamics

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/KCEvers/affectOU/graph/badge.svg)](https://app.codecov.io/gh/KCEvers/affectOU)
[![R-CMD-check](https://github.com/KCEvers/affectOU/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KCEvers/affectOU/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The aim of `affectOU` is to provide tools for using the
Ornstein-Uhlenbeck (OU) process in R with a chief focus on capturing
affect dynamics, or how feelings change over time. The OU process is
widely used for this purpose (e.g., Guthier et al., 2020; Kuppens et
al., 2010a; Voelkle & Oud, 2013), formulating three core psychological
mechanisms governed by three main sets of parameters (Uhlenbeck &
Ornstein, 1930; Oravecz et al., 2011):

- **Drift matrix $\Theta$:** Governs the rate at which affect returns to
  baseline, therefore capturing a person’s *emotion regulation capacity*
  or *emotional inertia*;
- **Attractor $\mu$:** Represents the long-term average level of affect
  that the system tends to return to, therefore capturing a person’s
  *baseline mood* or *emotional setpoint*;
- **Diffusion matrix $\Gamma$:** Controls the magnitude of short-term
  fluctuations or noise in the affective system, representing an
  individual’s sensitivity to environmental perturbations.

Rather than being a full-fledged package itself, `affectOU` rather
serves as a demonstration of how on can package computational models, as
detailed in Evers & Vanhasbroeck (2026).

## Installation

`affectOU` can be installed from
[GitHub](https://github.com/KCEvers/affectOU) with:

``` r
# install.packages("pak")
pak::pak("KCEvers/affectOU")
```

Once installed, one can use the package through calling the `library()`
function.

``` r
library(affectOU)
```

## Quick Start

The chief functionalities of the package revolve around a user-specified
OU model of a particular dimensionality. To create such a model, one
should use the
[`affectOU()`](https://kcevers.github.io/affectOU/reference/affectOU.html)
function and specify its dimensionality and/or parameters (see [Get
Started](https://kcevers.github.io/affectOU/articles/affectOU.html)).
For example, the following call creates a two-dimensional OU model with
default parameter values:

``` r
# Create two-dimensional OU
model <- affectOU(ndim = 2)
model
#> 
#> ── 2D Ornstein-Uhlenbeck Model ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> dX(t) = Θ(μ − X(t))dt + Γ dW(t)
#> 
#> μ = [0.000, 0.000]
#> 
#> Θ:
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> Γ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> Σ = ΓΓᵀ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```

To simulate data from this model, one should use the
[`simulate()`]((https://kcevers.github.io/affectOU/reference/simulate.affectOU.html))
function:

``` r
# Simulate timeseries
sim <- simulate(model)
```

To visualize this simulation, one uses the [`plot()`]() function
specifying the type of plot one wishes to create. For example, for a
time-series plot, one specifies `type = "time"`:

``` r
plot(sim, type = "time")
```

<img src="man/figures/README-unnamed-chunk-2-1.svg" alt="Visualization of a time-series plot, showing the results of the simulation. In one time-series plot, it shows the simulated values for both dimensions of the OU process together with a horizontal line denoting the baseline affective state $\mu$." width="100%" />

## Explore the OU Process

For more detailed examples and visual demonstrations of the model’s
characteristics, see the vignettes:

- [Get
  Started](https://kcevers.github.io/affectOU/articles/affectOU.html)
- [Affect Characteristics Implied by the OU
  Process](https://kcevers.github.io/affectOU/articles/characteristics.html)

Theoretical properties of the OU process are furthermore explained in
the documentation for the functions
[`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.html),
[`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.html),
and
[`relaxation()`](https://kcevers.github.io/affectOU/reference/relaxation.affectOU.html).
For fitting unidimensional OU models to data, see
[`fit()`](https://kcevers.github.io/affectOU/reference/fit.affectOU.html).

## References

Evers, K. C. & Vanhasbroeck, N. (2026). *Sharing computational models as
reproducible and user-friendly packages: A tutorial in R.* Manuscript in
preparation.

Guthier, C., Dormann, C., & Voelkle, M. C. (2020). Reciprocal effects
between job stressors and burnout: A continuous-time meta-analysis of
longitudinal studies. *Psychological Bulletin, 146*, 1146-1173. doi:
[10.1037/bul0000304](https://doi.org/10.1037/bul0000304)

Kuppens, P., Oravecz, Z., & Tuerlinckx, F. (2010). Feelings change:
Accounting for individual differences in the temporal dynamics of
affect. *Journal of Personality and Social Psychology, 99*, 1042-1060.
doi: [10.1037/a0020962](https://doi.org/10.1037/a0020962)

Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011). A hierarchical
latent stochastic differential equation model for affective dynamics.
*Psychological Methods, 16*, 468-490. doi:
[10.1037/a0024375](https://doi.org/10.1037/a0024375)

Uhlenbeck, G. E. & Ornstein, L. S. (1930). On the theory of Brownian
motion. *Physical Review, 46*, 823-841.

Voelkle, M. C. & Oud, J. H. L. (2013). Continuous time modelling with
individually varying time intervals for oscillating and non-oscillating
processes. *British Journal of Mathematical and Statistical Psychology,
66*, 103-126. doi:
[10.1111/j.2044-8317.2012.02043.x](https://doi.org/10.1111/j.2044-8317.2012.02043.x)
