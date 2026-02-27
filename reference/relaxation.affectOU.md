# Compute relaxation time for Ornstein-Uhlenbeck model

The relaxation time \\\tau_k = 1/\mathrm{Re}(\lambda_k)\\ is the
characteristic time scale of eigenmode \\k\\ of the OU process. It
measures how quickly deviations from equilibrium decay. For oscillatory
modes (complex eigenvalues) the oscillation period \\T_k = 2\pi /
\|\mathrm{Im}(\lambda_k)\|\\ is also reported alongside.

## Usage

``` r
# S3 method for class 'affectOU'
relaxation(object, ...)
```

## Arguments

- object:

  An `affectOU` model object.

- ...:

  Additional arguments (unused).

## Value

A data frame of class `relaxation_affectOU` with one row per eigenmode.
Complex conjugate pairs are collapsed to a single oscillatory mode.
Columns:

- `mode`:

  Integer mode index

- `relaxation_time`:

  Relaxation time \\\tau_k = 1/\mathrm{Re}(\lambda_k)\\, or `NA` if
  \\\mathrm{Re}(\lambda_k) \le 0\\ (unstable mode)

- `half_life`:

  Half-life \\t\_{1/2} = \ln 2 \cdot \tau_k\\, or `NA`

- `oscillation_period`:

  Period \\T_k = 2\pi/\|\mathrm{Im}(\lambda_k)\|\\ for oscillatory
  modes, `NA` for real eigenvalues

- `eigenvalue_re`:

  Real part of eigenvalue \\\mathrm{Re}(\lambda_k)\\

- `eigenvalue_im`:

  Imaginary part magnitude \\\|\mathrm{Im}(\lambda_k)\|\\ (0 for real
  eigenvalues)

- `is_oscillatory`:

  `TRUE` for complex eigenvalue pairs

Attributes:

- `tau_max`:

  Slowest mode: \\1/\min_k \mathrm{Re}(\lambda_k)\\ over stable modes

- `tau_min`:

  Fastest mode: \\1/\max_k \mathrm{Re}(\lambda_k)\\ over stable modes

- `half_life_max`:

  \\\ln 2 \cdot \tau\_{\max}\\

- `half_life_min`:

  \\\ln 2 \cdot \tau\_{\min}\\

- `ndim`:

  Dimensionality of the process

## Exponential decay

For a stable 1D OU process starting from \\x_0\\, the expected
trajectory is: \$\$E\[X(t)\] = \mu + (x_0 - \mu)\\e^{-\theta t}\$\$ The
deviation from baseline decays exponentially. The relaxation time (also
called *decorrelation time*) is the e-folding time: \$\$\tau =
\frac{1}{\theta}\$\$ At lag \\\tau\\, the autocorrelation function
\\\mathrm{ACF}(\tau) = e^{-\theta\tau} = 1/e \approx 0.368\\, which is
exactly the decorrelation time. The half-life is the lag at which 50\\
\$\$t\_{1/2} = \ln 2 \cdot \tau = \frac{\ln 2}{\theta} \approx
\frac{0.693}{\theta}\$\$

## Eigenvalue-based relaxation time (multivariate)

For a multivariate OU process, the expected trajectory from perturbation
\\x_0\\ is: \$\$E\[X(t) \mid X(0) = x_0\] = \mu + e^{-\Theta t}(x_0 -
\mu)\$\$ where \\e^{-\Theta t}\\ is the matrix exponential. Via the
eigendecomposition \\\Theta = P \Lambda P^{-1}\\: \$\$e^{-\Theta t} =
P\\\mathrm{diag}(e^{-\lambda_1 t},\ldots, e^{-\lambda_n t})\\P^{-1}\$\$

Each eigenmode \\k\\ decays with time constant \\\tau_k =
1/\mathrm{Re}(\lambda_k)\\ (defined only when \\\mathrm{Re}(\lambda_k)
\> 0\\). Key summary timescales:

- **Slowest mode**: \\\tau\_{\max} = 1/\min_k \mathrm{Re}(\lambda_k)\\ —
  dominates long-lag behaviour and is often the most psychologically
  interpretable timescale (emotional inertia).

- **Fastest mode**: \\\tau\_{\min} = 1/\max_k \mathrm{Re}(\lambda_k)\\ —
  governs the fastest transient dynamics.

## Decorrelation time

The relaxation time \\\tau = 1/\mathrm{Re}(\lambda)\\ is also the
*decorrelation time*: the e-folding lag of the autocorrelation function.
For the 1D case, \\\mathrm{ACF}(\tau) = e^{-\theta\tau} = 1/e\\ exactly.
For coupled multivariate systems, the long-lag ACF of every dimension is
dominated by the slowest eigenmode, so \\\tau\_{\max} =
1/\mathrm{Re}(\lambda\_{\min})\\ is the system-level decorrelation time.

## Oscillatory modes (stable spiral)

When \\\Theta\\ has complex eigenvalues \\\lambda_k = a \pm bi\\ (\\a \>
0\\), the corresponding mode produces *damped oscillations*: the
perturbation decays with envelope time constant \\\tau_k = 1/a\\ while
oscillating with period \\T_k = 2\pi/b\\. Conjugate pairs share the same
\\\tau_k\\ and \\T_k\\ and are reported as a single oscillatory mode.

## See also

[`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)
for stability assessment,
[`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md)
for the equilibrium distribution,
[`summary()`](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
for the full model summary.

## Examples

``` r
# 1D stable
model <- affectOU(theta = 0.5)
relaxation(model)
#> 
#> ── Relaxation time of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Relaxation time (τ): 2 time units
#> Half-life (t₁/₂): 1.386 time units
#> 
#> Time for perturbation to decay to:
#>      50%     37%     14%      5%      1%
#>    1.386   2.000   4.000   6.000  10.000

# 1D random walk (no relaxation)
model_rw <- affectOU(theta = 0)
relaxation(model_rw)
#> 
#> ── Relaxation time of 1D Ornstein-Uhlenbeck Model ──
#> 
#> All modes are unstable (no finite relaxation time).

# 2D diagonal: two independent modes
model_diag <- affectOU(theta = diag(c(0.5, 0.2)), mu = 0, gamma = 1)
relaxation(model_diag)
#> 
#> ── Relaxation time of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Slowest mode: τ_max = 5 (t₁/₂ = 3.466) time units
#> Fastest mode: τ_min = 2 (t₁/₂ = 1.386) time units
#> 
#> Mode 1: τ = 5, t₁/₂ = 3.466 (λ = 0.2)
#> Mode 2: τ = 2, t₁/₂ = 1.386 (λ = 0.5)
#> 
#> Time for perturbation to decay to:
#>               50%     37%     14%      5%      1%
#>   Mode 1:   3.466   5.000  10.000  15.000  25.000
#>   Mode 2:   1.386   2.000   4.000   6.000  10.000

# 2D stable spiral: one oscillatory mode
theta_sp <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
model_sp <- affectOU(theta = theta_sp, mu = 0, gamma = 1)
relaxation(model_sp)
#> 
#> ── Relaxation time of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Relaxation time (τ): 2 time units
#> Half-life (t₁/₂): 1.386 time units
#> 
#> Time for perturbation to decay to (envelope for oscillatory modes):
#>      50%     37%     14%      5%      1%
#>    1.386   2.000   4.000   6.000  10.000

# 2D mixed stability: one stable, one unstable mode
model_mixed <- suppressWarnings(affectOU(theta = diag(c(0.5, -0.3))))
relaxation(model_mixed)
#> 
#> ── Relaxation time of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Slowest mode: τ_max = 2 (t₁/₂ = 1.386) time units
#> Fastest mode: τ_min = 2 (t₁/₂ = 1.386) time units
#> 
#> Mode 1 (unstable): τ = NA (λ = -0.3)
#> Mode 2: τ = 2, t₁/₂ = 1.386 (λ = 0.5)
#> 
#> Time for perturbation to decay to:
#>               50%     37%     14%      5%      1%
#>   Mode 2:   1.386   2.000   4.000   6.000  10.000
```
