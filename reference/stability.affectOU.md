# Assess stability of an Ornstein-Uhlenbeck model

Classify the dynamics and stability of an Ornstein-Uhlenbeck model based
on eigenvalue analysis of the drift matrix \\\Theta\\.

## Usage

``` r
# S3 method for class 'affectOU'
stability(object, tol = 1e-10, ...)
```

## Arguments

- object:

  An `affectOU` model object.

- tol:

  Tolerance for comparing eigenvalues to zero.

- ...:

  Additional arguments (unused).

## Value

A list of class `stability_affectOU` containing:

- is_stable:

  Logical, `TRUE` if the system is stable (all eigenvalues have positive
  real parts)

- dynamics:

  Character string describing the dynamics type using dynamical systems
  terminology. For 1D: `"stable node"`, `"random walk"`, or
  `"unstable node"`. For multivariate: `"stable node"`,
  `"stable spiral"`, `"unstable node"`, `"unstable spiral"`,
  `"marginally stable"`, `"saddle point"`, or `"saddle spiral"`.

- per_dimension:

  Character vector with per-dimension dynamics classification (`NULL`
  for 1D)

- eigenvalues:

  Eigenvalues of theta

- ndim:

  Dimensionality of the process

## Details

Note that although typically, a system is stable if all its eigenvalues
have negative real parts, due to the parameterization of the OU process,
stability is determined by *positive* real parts:

\$\$dX(t) = \Theta(\mu - X(t))dt + \Gamma dW(t)\$\$

As \\\Theta\\ is multiplied by the negative of the state variable,
positive eigenvalues indicate that deviations from the attractor decay
over time, leading to a stable system. Conversely, negative eigenvalues
indicate that deviations grow over time, leading to an unstable system.

## Stability via eigenvalues

A stationary distribution exists if and only if all eigenvalues of
\\\Theta\\ have positive real parts. This is a system-level property,
not a per-dimension property:

1.  Positive diagonal elements (\\\Theta\_{ii} \> 0\\) alone do not
    guarantee stability. Strong off-diagonal coupling can push
    eigenvalues toward zero or negative real parts, destabilising the
    system.

2.  Conversely, coupling from other dimensions can stabilise a dimension
    whose diagonal element is small or even zero. A dimension that would
    be non-stationary in isolation may become stationary when embedded
    in a coupled system.

## Oscillatory dynamics

Complex eigenvalues (arising from asymmetric coupling) produce
oscillatory dynamics, where dimensions cycle around each other. When all
eigenvalues have positive real parts, these oscillations are damped and
the system still converges to \\\mu\\. When any eigenvalue has a
non-positive real part, the oscillations grow and the system diverges.

## Terminology

System-level dynamics are classified using dynamical systems terms: a
node converges or diverges without oscillation, a spiral exhibits damped
or growing oscillations, and a saddle point has mixed
convergence/divergence across directions. Per-dimension dynamics
describe the behavior of each dimension in isolation based on its
diagonal element in \\\Theta\\ ("stable node" if positive, "random walk"
if zero, "unstable node" if negative). A dimension classified as "stable
node" in isolation may still belong to an unstable coupled system due to
cross-regulation.

## See also

[stationary()](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md)
for the equilibrium distribution,
[relaxation()](https://kcevers.github.io/affectOU/reference/relaxation.affectOU.md)
for perturbation persistence,
[summary()](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
for the full model summary

## Examples

``` r
# 1D stable node
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
stability(model)
#> 
#> ── Stability analysis of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Stable (node). Deviations from the attractor decay exponentially.

# 1D random walk (not stable)
model_rw <- affectOU(theta = 0)
stability(model_rw)
#> 
#> ── Stability analysis of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Not stable (random walk). No attractor; the process drifts freely.

# Positive diagonals with oscillatory coupling: stable
theta_osc <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
eigen(theta_osc)$values # complex with positive real part
#> [1] 0.5+0.4i 0.5-0.4i
stability(affectOU(theta = theta_osc, mu = 0, gamma = 1))
#> 
#> ── Stability analysis of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Stable (spiral). The system spirals toward the attractor with damped
#> oscillations.

# Strong coupling still stable if real parts stay positive
theta_strong <- matrix(c(0.5, -1.5, 1.5, 0.5), nrow = 2)
stability(affectOU(theta = theta_strong, mu = 0, gamma = 1))
#> 
#> ── Stability analysis of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Stable (spiral). The system spirals toward the attractor with damped
#> oscillations.

# All diagonals positive, but coupling destabilises the system
theta_destab <- matrix(c(0.5, 1.0, 1.0, 0.5), nrow = 2)
stability(affectOU(theta = theta_destab, mu = 0, gamma = 1))
#> 
#> ── Stability analysis of 2D Ornstein-Uhlenbeck Model ──
#> 
#> Not stable (saddle point). Some directions converge while others diverge.

# One negative diagonal element makes the system non-stationary
theta_unstable <- matrix(c(0.5, 0, 0, -0.3), nrow = 2)
stability(affectOU(theta = theta_unstable, mu = 0, gamma = 1))$is_stable
#> [1] FALSE
```
