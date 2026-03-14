# print.affectOU snapshot (1D) [plain]

    Code
      print(model)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = θ(μ − X(t))dt + γ dW(t)
      
      θ = 0.500, μ = 0.000, γ = 1.000, σ = γ² = 1.000

---

    Code
      print(model, digits = 5)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = θ(μ − X(t))dt + γ dW(t)
      
      θ = 0.50000, μ = 0.00000, γ = 1.00000, σ = γ² = 1.00000

# print.affectOU snapshot (1D) [ansi]

    Code
      print(model)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[3mX[23m(t) = [3mθ[23m([3mμ[23m − [3mX[23m(t))dt + [3mγ[23m d[3mW[23m(t)
      
      [3mθ[23m = 0.500, [3mμ[23m = 0.000, [3mγ[23m = 1.000, [3mσ[23m = [3mγ[23m² = 1.000

---

    Code
      print(model, digits = 5)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[3mX[23m(t) = [3mθ[23m([3mμ[23m − [3mX[23m(t))dt + [3mγ[23m d[3mW[23m(t)
      
      [3mθ[23m = 0.50000, [3mμ[23m = 0.00000, [3mγ[23m = 1.00000, [3mσ[23m = [3mγ[23m² = 1.00000

# print.affectOU snapshot (2D) [plain]

    Code
      print(model)
    Message
      
      -- 2D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = Θ(μ − X(t))dt + Γ dW(t)
      
      μ = [1.000, -1.000]
      
      Θ:
           [,1] [,2]
      [1,]  0.5  0.0
      [2,]  0.0  0.3
      
      Γ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1
      
      Σ = ΓΓᵀ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1

---

    Code
      print(model, digits = 5)
    Message
      
      -- 2D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = Θ(μ − X(t))dt + Γ dW(t)
      
      μ = [1.00000, -1.00000]
      
      Θ:
           [,1] [,2]
      [1,]  0.5  0.0
      [2,]  0.0  0.3
      
      Γ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1
      
      Σ = ΓΓᵀ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1

# print.affectOU snapshot (2D) [ansi]

    Code
      print(model)
    Message
      
      [36m--[39m [1m2D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[1mX[22m(t) = [1mΘ[22m([1mμ[22m − [1mX[22m(t))dt + [1mΓ[22m d[1mW[22m(t)
      
      [1mμ[22m = [1.000, -1.000]
      
      [1mΘ[22m:
           [,1] [,2]
      [1,]  0.5  0.0
      [2,]  0.0  0.3
      
      [1mΓ[22m:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1
      
      [1mΣ[22m = [1mΓΓ[22mᵀ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1

---

    Code
      print(model, digits = 5)
    Message
      
      [36m--[39m [1m2D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[1mX[22m(t) = [1mΘ[22m([1mμ[22m − [1mX[22m(t))dt + [1mΓ[22m d[1mW[22m(t)
      
      [1mμ[22m = [1.00000, -1.00000]
      
      [1mΘ[22m:
           [,1] [,2]
      [1,]  0.5  0.0
      [2,]  0.0  0.3
      
      [1mΓ[22m:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1
      
      [1mΣ[22m = [1mΓΓ[22mᵀ:
           [,1] [,2]
      [1,]    1    0
      [2,]    0    1

# print.affectOU snapshot (21D) [plain]

    Code
      print(model)
    Message
      
      -- 21D Ornstein-Uhlenbeck Model ------------------------------------------------
      dX(t) = Θ(μ − X(t))dt + Γ dW(t)
      
      i High-dimensional model. Parameters are not shown, but can be accessed with `coef()`.

# print.affectOU snapshot (21D) [ansi]

    Code
      print(model)
    Message
      
      [36m--[39m [1m21D Ornstein-Uhlenbeck Model[22m [36m------------------------------------------------[39m
      d[1mX[22m(t) = [1mΘ[22m([1mμ[22m − [1mX[22m(t))dt + [1mΓ[22m d[1mW[22m(t)
      
      [36mi[39m High-dimensional model. Parameters are not shown, but can be accessed with `coef()`.

# print.affectOU digits snapshot [plain]

    Code
      print(model, digits = 1)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = θ(μ − X(t))dt + γ dW(t)
      
      θ = 0.1, μ = 0.0, γ = 1.0, σ = γ² = 1.0

---

    Code
      print(model, digits = 5)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      dX(t) = θ(μ − X(t))dt + γ dW(t)
      
      θ = 0.12346, μ = 0.00000, γ = 1.00000, σ = γ² = 1.00000

# print.affectOU digits snapshot [ansi]

    Code
      print(model, digits = 1)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[3mX[23m(t) = [3mθ[23m([3mμ[23m − [3mX[23m(t))dt + [3mγ[23m d[3mW[23m(t)
      
      [3mθ[23m = 0.1, [3mμ[23m = 0.0, [3mγ[23m = 1.0, [3mσ[23m = [3mγ[23m² = 1.0

---

    Code
      print(model, digits = 5)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      d[3mX[23m(t) = [3mθ[23m([3mμ[23m − [3mX[23m(t))dt + [3mγ[23m d[3mW[23m(t)
      
      [3mθ[23m = 0.12346, [3mμ[23m = 0.00000, [3mγ[23m = 1.00000, [3mσ[23m = [3mγ[23m² = 1.00000

# print.summary_affectOU snapshot (1D) [plain]

    Code
      print(s)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      
      -- Dynamics --
      
      Stable (node)
      
      -- Stationary distribution --
      
      Mean: 0
      SD: 1

# print.summary_affectOU snapshot (1D) [ansi]

    Code
      print(s)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      
      -- [1m[1mDynamics[1m[22m --
      
      Stable (node)
      
      -- [1m[1mStationary distribution[1m[22m --
      
      Mean: 0
      SD: 1

# print.summary_affectOU snapshot (2D) [plain]

    Code
      print(s)
    Message
      
      -- 2D Ornstein-Uhlenbeck Model -------------------------------------------------
      
      -- Dynamics --
      
      Stable (node)
      
      -- Stationary distribution --
      
      Mean: [0, 0]
      SD: [1, 1.291]
      
      -- Structure --
      
      Coupling: none
      Noise: independent

# print.summary_affectOU snapshot (2D) [ansi]

    Code
      print(s)
    Message
      
      [36m--[39m [1m2D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      
      -- [1m[1mDynamics[1m[22m --
      
      Stable (node)
      
      -- [1m[1mStationary distribution[1m[22m --
      
      Mean: [0, 0]
      SD: [1, 1.291]
      
      -- [1m[1mStructure[1m[22m --
      
      Coupling: none
      Noise: independent

# print.summary_affectOU snapshot (coupled) [plain]

    Code
      print(s)
    Message
      
      -- 2D Ornstein-Uhlenbeck Model -------------------------------------------------
      
      -- Dynamics --
      
      Stable (node)
      
      -- Stationary distribution --
      
      Mean: [0, 0]
      SD: [1.035, 1.336]
      
      -- Structure --
      
      Coupling: Dim 1 → Dim 2 (+), Dim 2 → Dim 1 (+)
      Noise: independent

# print.summary_affectOU snapshot (coupled) [ansi]

    Code
      print(s)
    Message
      
      [36m--[39m [1m2D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      
      -- [1m[1mDynamics[1m[22m --
      
      Stable (node)
      
      -- [1m[1mStationary distribution[1m[22m --
      
      Mean: [0, 0]
      SD: [1.035, 1.336]
      
      -- [1m[1mStructure[1m[22m --
      
      Coupling: Dim 1 → Dim 2 (+), Dim 2 → Dim 1 (+)
      Noise: independent

# print.summary_affectOU snapshot (unstable) [plain]

    Code
      print(s)
    Message
      
      -- 1D Ornstein-Uhlenbeck Model -------------------------------------------------
      
      -- Dynamics --
      
      Not stable (node)
      
      -- Stationary distribution --
      
      Does not exist (system is not stable).

# print.summary_affectOU snapshot (unstable) [ansi]

    Code
      print(s)
    Message
      
      [36m--[39m [1m1D Ornstein-Uhlenbeck Model[22m [36m-------------------------------------------------[39m
      
      -- [1m[1mDynamics[1m[22m --
      
      Not stable (node)
      
      -- [1m[1mStationary distribution[1m[22m --
      
      Does not exist (system is not stable).

