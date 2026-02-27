# relaxation print: 1D stable [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 1D Ornstein-Uhlenbeck Model --
      
      Relaxation time (τ): 2 time units
      Half-life (t₁/₂): 1.386 time units
      
      Time for perturbation to decay to:
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

# relaxation print: 1D stable [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Relaxation time (τ): 2 time units
      Half-life (t₁/₂): 1.386 time units
      
      Time for perturbation to decay to:
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

# relaxation print: 1D unstable [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 1D Ornstein-Uhlenbeck Model --
      
      All modes are unstable (no finite relaxation time).

# relaxation print: 1D unstable [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      All modes are unstable (no finite relaxation time).

# relaxation print: 2D diagonal [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 2D Ornstein-Uhlenbeck Model --
      
      Slowest mode: τ_max = 5 (t₁/₂ = 3.466) time units
      Fastest mode: τ_min = 2 (t₁/₂ = 1.386) time units
      
      Mode 1: τ = 5, t₁/₂ = 3.466 (λ = 0.2)
      Mode 2: τ = 2, t₁/₂ = 1.386 (λ = 0.5)
      
      Time for perturbation to decay to:
                    50%     37%     14%      5%      1%
        Mode 1:   3.466   5.000  10.000  15.000  25.000
        Mode 2:   1.386   2.000   4.000   6.000  10.000

# relaxation print: 2D diagonal [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Slowest mode: τ_max = 5 (t₁/₂ = 3.466) time units
      Fastest mode: τ_min = 2 (t₁/₂ = 1.386) time units
      
      Mode 1: τ = 5, t₁/₂ = 3.466 (λ = 0.2)
      Mode 2: τ = 2, t₁/₂ = 1.386 (λ = 0.5)
      
      Time for perturbation to decay to:
                    50%     37%     14%      5%      1%
        Mode 1:   3.466   5.000  10.000  15.000  25.000
        Mode 2:   1.386   2.000   4.000   6.000  10.000

# relaxation print: 2D stable spiral [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 2D Ornstein-Uhlenbeck Model --
      
      Relaxation time (τ): 2 time units
      Half-life (t₁/₂): 1.386 time units
      
      Time for perturbation to decay to (envelope for oscillatory modes):
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

# relaxation print: 2D stable spiral [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Relaxation time (τ): 2 time units
      Half-life (t₁/₂): 1.386 time units
      
      Time for perturbation to decay to (envelope for oscillatory modes):
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

