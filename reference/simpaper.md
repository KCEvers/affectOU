# Simulation in Evers & Vanhasbroeck (2026)

A two-dimensional Ornstein-Uhlenbeck process simulated with `affectOU`.

## Usage

``` r
simpaper
```

## Format

An object of class `simulate_affectOU` containing simulated trajectories
from the specified Ornstein-Uhlenbeck process.

## Details

The process is defined by the stochastic differential equation:
\$\$d\mathbf{X}\_t = \mathbf{\Theta}(\mathbf{\mu} - \mathbf{X}\_t)dt +
\mathbf{\Gamma}^{1/2}d\mathbf{W}\_t\$\$

where:

- \\\mathbf{\Theta} = \begin{pmatrix} 0.7 & 0 \\ 0 & 0.3 \end{pmatrix}\\
  is the drift matrix

- \\\mathbf{\mu} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}\\ is the
  equilibrium mean vector

- \\\mathbf{\Gamma} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\\ is
  the diffusion matrix

- \\\mathbf{W}\_t\\ is a two-dimensional Wiener process
