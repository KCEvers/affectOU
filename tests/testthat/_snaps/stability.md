# stability print shows 1D stable node [plain]

    Code
      print(stability(model))
    Message
      
      -- Stability analysis of 1D Ornstein-Uhlenbeck Model --
      
      Stable (node). Deviations from the attractor decay exponentially.
      
      Eigenvalues (all real):
      * λ1: 0.5

# stability print shows 1D stable node [ansi]

    Code
      print(stability(model))
    Message
      
      -- [1m[1mStability analysis of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Stable (node). Deviations from the attractor decay exponentially.
      
      Eigenvalues (all real):
      * λ1: 0.5

# stability print shows 1D non-stationary [plain]

    Code
      print(stability(model))
    Message
      
      -- Stability analysis of 1D Ornstein-Uhlenbeck Model --
      
      Not stable (node). Deviations from the attractor grow exponentially.
      
      Eigenvalues (all real):
      * λ1: -0.5

# stability print shows 1D non-stationary [ansi]

    Code
      print(stability(model))
    Message
      
      -- [1m[1mStability analysis of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Not stable (node). Deviations from the attractor grow exponentially.
      
      Eigenvalues (all real):
      * λ1: -0.5

# stability print shows 2D stable [plain]

    Code
      print(stability(model))
    Message
      
      -- Stability analysis of 2D Ornstein-Uhlenbeck Model --
      
      Stable (node). Deviations from the attractor decay exponentially.
      
      Eigenvalues (all real):
      * λ1: 0.5
      * λ2: 0.3

# stability print shows 2D stable [ansi]

    Code
      print(stability(model))
    Message
      
      -- [1m[1mStability analysis of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Stable (node). Deviations from the attractor decay exponentially.
      
      Eigenvalues (all real):
      * λ1: 0.5
      * λ2: 0.3

# stability print shows 2D oscillatory [plain]

    Code
      print(stability(model))
    Message
      
      -- Stability analysis of 2D Ornstein-Uhlenbeck Model --
      
      Stable (spiral). The system spirals toward the attractor with damped
      oscillations.
      
      Eigenvalues (complex):
      * λ1: 0.5 + 0.4i
      * λ2: 0.5 - 0.4i

# stability print shows 2D oscillatory [ansi]

    Code
      print(stability(model))
    Message
      
      -- [1m[1mStability analysis of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Stable (spiral). The system spirals toward the attractor with damped
      oscillations.
      
      Eigenvalues (complex):
      * λ1: 0.5 + 0.4i
      * λ2: 0.5 - 0.4i

