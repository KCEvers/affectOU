withr::with_seed(
  seed = 123,
  code = {
    model <- affectOU(ndim = 2, theta = 1, mu = 0, gamma = 1)
    simpaper <- simulate(model)
    usethis::use_data(simpaper, overwrite = TRUE)
  },
  .rng_kind = "Mersenne-Twister",
  .rng_normal_kind = "Inversion",
  .rng_sample_kind = "Rejection"
)
