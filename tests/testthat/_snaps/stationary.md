# stationary print shows 1D model [plain]

    Code
      print(stationary(model))
    Message
      
      -- Stationary distribution of 1D Ornstein-Uhlenbeck Model --
      
      Mean: 0
      SD: 1
      95% interval: [-2, 2]

# stationary print shows 1D model [ansi]

    Code
      print(stationary(model))
    Message
      
      -- [1m[1mStationary distribution of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Mean: 0
      SD: 1
      95% interval: [-2, 2]

# stationary print shows non-stable model [plain]

    Code
      print(stationary(model))
    Message
      
      -- Stationary distribution of 1D Ornstein-Uhlenbeck Model --
      
      Does not exist (system is not stable).

# stationary print shows non-stable model [ansi]

    Code
      print(stationary(model))
    Message
      
      -- [1m[1mStationary distribution of [1m1[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Does not exist (system is not stable).

# stationary print shows 2D model [plain]

    Code
      print(stationary(model))
    Message
      
      -- Stationary distribution of 2D Ornstein-Uhlenbeck Model --
      
      Mean: [0, 0]
      SD: [1.035, 1.336]
      
      Stationary correlations:
      * Dims 1 & 2: -0.258

# stationary print shows 2D model [ansi]

    Code
      print(stationary(model))
    Message
      
      -- [1m[1mStationary distribution of [1m2[1mD Ornstein-Uhlenbeck Model[1m[22m --
      
      Mean: [0, 0]
      SD: [1.035, 1.336]
      
      Stationary correlations:
      * Dims 1 & 2: -0.258

