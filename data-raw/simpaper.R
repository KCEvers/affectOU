withr::with_seed(
  seed = 123,
  code = {
    theta_2d <- matrix(c(0.7, 0, 0, 0.3), nrow = 2)
    model <- affectOU(theta = theta_2d, mu = 0, gamma = 1)
    simpaper <- simulate(model, seed = 123)
    usethis::use_data(simpaper, overwrite = TRUE)
  },
  .rng_kind = "Mersenne-Twister",
  .rng_normal_kind = "Inversion",
  .rng_sample_kind = "Rejection"
)
