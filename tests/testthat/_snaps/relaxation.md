# relaxation print shows heading and bullets for 1D model [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 1D Ornstein-Uhlenbeck Model --
      
      Half-life (t₁/₂): 1.386 time units
      Relaxation time (τ): 2 time units
      
      Time for perturbation to decay to:
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

# relaxation print shows heading and bullets for 1D model [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Half-life (t₁/₂): 1.386 time units
      Relaxation time (τ): 2 time units
      
      Time for perturbation to decay to:
           50%     37%     14%      5%      1%
         1.386   2.000   4.000   6.000  10.000

# relaxation print handles NA relaxation time gracefully [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 1D Ornstein-Uhlenbeck Model --
      
      Relaxation time: undefined (not stable)

# relaxation print handles NA relaxation time gracefully [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Relaxation time: undefined (not stable)

# relaxation print shows multi-dim result [plain]

    Code
      print(relaxation(model))
    Message
      
      -- Relaxation time of 2D Ornstein-Uhlenbeck Model --
      
      Dim. 1: t₁/₂ = 1.386, τ = 2
      Dim. 2: t₁/₂ = 3.466, τ = 5
      
      Time for perturbation to decay to:
                    50%     37%     14%      5%      1%
        Dim. 1:   1.386   2.000   4.000   6.000  10.000
        Dim. 2:   3.466   5.000  10.000  15.000  25.000

# relaxation print shows multi-dim result [ansi]

    Code
      print(relaxation(model))
    Message
      
      -- [1m[1mRelaxation time of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Dim. 1: t₁/₂ = 1.386, τ = 2
      Dim. 2: t₁/₂ = 3.466, τ = 5
      
      Time for perturbation to decay to:
                    50%     37%     14%      5%      1%
        Dim. 1:   1.386   2.000   4.000   6.000  10.000
        Dim. 2:   3.466   5.000  10.000  15.000  25.000

