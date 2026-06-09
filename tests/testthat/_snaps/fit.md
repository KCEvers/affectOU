# print.fit_affectOU snapshot [plain]

    Code
      print(fit_obj)
    Message
      
      -- Fitted 1D Ornstein-Uhlenbeck Model ------------------------------------------
      101 data points (dt ≈ 0.100)
      θ = 0.740, μ = 0.249, γ = 1.071
      Log-likelihood: -29.935
      RMSE: 0.325

# print.fit_affectOU snapshot [ansi]

    Code
      print(fit_obj)
    Message
      
      [36m--[39m [1mFitted 1D Ornstein-Uhlenbeck Model[22m [36m------------------------------------------[39m
      101 data points (dt ≈ 0.100)
      [3mθ[23m = 0.740, [3mμ[23m = 0.249, [3mγ[23m = 1.071
      Log-likelihood: -29.935
      RMSE: 0.325

# print.fit_affectOU digits snapshot [plain]

    Code
      print(fit_obj, digits = 1)
    Message
      
      -- Fitted 1D Ornstein-Uhlenbeck Model ------------------------------------------
      101 data points (dt ≈ 0.1)
      θ = 0.7, μ = 0.2, γ = 1.1
      Log-likelihood: -29.9
      RMSE: 0.3

---

    Code
      print(fit_obj, digits = 5)
    Message
      
      -- Fitted 1D Ornstein-Uhlenbeck Model ------------------------------------------
      101 data points (dt ≈ 0.10000)
      θ = 0.74024, μ = 0.24944, γ = 1.07065
      Log-likelihood: -29.93467
      RMSE: 0.32479

# print.fit_affectOU digits snapshot [ansi]

    Code
      print(fit_obj, digits = 1)
    Message
      
      [36m--[39m [1mFitted 1D Ornstein-Uhlenbeck Model[22m [36m------------------------------------------[39m
      101 data points (dt ≈ 0.1)
      [3mθ[23m = 0.7, [3mμ[23m = 0.2, [3mγ[23m = 1.1
      Log-likelihood: -29.9
      RMSE: 0.3

---

    Code
      print(fit_obj, digits = 5)
    Message
      
      [36m--[39m [1mFitted 1D Ornstein-Uhlenbeck Model[22m [36m------------------------------------------[39m
      101 data points (dt ≈ 0.10000)
      [3mθ[23m = 0.74024, [3mμ[23m = 0.24944, [3mγ[23m = 1.07065
      Log-likelihood: -29.93467
      RMSE: 0.32479

# print.summary_fit_affectOU snapshot [plain]

    Code
      print(s)
    Message
      
      -- Fitted Ornstein-Uhlenbeck Model Summary -------------------------------------
      Method: mle, 101 observations
      
      -- Coefficients --
      
            Estimate    SE          95% CI
      theta    0.740 0.382 [-0.009, 1.489]
      mu       0.249 0.462 [-0.656, 1.155]
      gamma    1.071 0.078  [0.917, 1.224]
      
      -- Goodness of fit --
      
      Log-likelihood: -29.935
      RMSE: 0.325

# print.summary_fit_affectOU snapshot [ansi]

    Code
      print(s)
    Message
      
      [36m--[39m [1mFitted Ornstein-Uhlenbeck Model Summary[22m [36m-------------------------------------[39m
      Method: mle, 101 observations
      
      -- [1m[1mCoefficients[1m[22m --
      
            Estimate    SE          95% CI
      theta    0.740 0.382 [-0.009, 1.489]
      mu       0.249 0.462 [-0.656, 1.155]
      gamma    1.071 0.078  [0.917, 1.224]
      
      -- [1m[1mGoodness of fit[1m[22m --
      
      Log-likelihood: -29.935
      RMSE: 0.325

# print.summary_fit_affectOU digits snapshot [plain]

    Code
      print(s, digits = 1)
    Message
      
      -- Fitted Ornstein-Uhlenbeck Model Summary -------------------------------------
      Method: mle, 101 observations
      
      -- Coefficients --
      
            Estimate  SE      95% CI
      theta      0.7 0.4 [-0.0, 1.5]
      mu         0.2 0.5 [-0.7, 1.2]
      gamma      1.1 0.1  [0.9, 1.2]
      
      -- Goodness of fit --
      
      Log-likelihood: -29.9
      RMSE: 0.3

---

    Code
      print(s, digits = 5)
    Message
      
      -- Fitted Ornstein-Uhlenbeck Model Summary -------------------------------------
      Method: mle, 101 observations
      
      -- Coefficients --
      
            Estimate      SE              95% CI
      theta  0.74024 0.38219 [-0.00884, 1.48933]
      mu     0.24944 0.46188 [-0.65581, 1.15470]
      gamma  1.07065 0.07829  [0.91720, 1.22410]
      
      -- Goodness of fit --
      
      Log-likelihood: -29.93467
      RMSE: 0.32479

# print.summary_fit_affectOU digits snapshot [ansi]

    Code
      print(s, digits = 1)
    Message
      
      [36m--[39m [1mFitted Ornstein-Uhlenbeck Model Summary[22m [36m-------------------------------------[39m
      Method: mle, 101 observations
      
      -- [1m[1mCoefficients[1m[22m --
      
            Estimate  SE      95% CI
      theta      0.7 0.4 [-0.0, 1.5]
      mu         0.2 0.5 [-0.7, 1.2]
      gamma      1.1 0.1  [0.9, 1.2]
      
      -- [1m[1mGoodness of fit[1m[22m --
      
      Log-likelihood: -29.9
      RMSE: 0.3

---

    Code
      print(s, digits = 5)
    Message
      
      [36m--[39m [1mFitted Ornstein-Uhlenbeck Model Summary[22m [36m-------------------------------------[39m
      Method: mle, 101 observations
      
      -- [1m[1mCoefficients[1m[22m --
      
            Estimate      SE              95% CI
      theta  0.74024 0.38219 [-0.00884, 1.48933]
      mu     0.24944 0.46188 [-0.65581, 1.15470]
      gamma  1.07065 0.07829  [0.91720, 1.22410]
      
      -- [1m[1mGoodness of fit[1m[22m --
      
      Log-likelihood: -29.93467
      RMSE: 0.32479

