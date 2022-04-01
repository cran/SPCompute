## Testing for Power Computation using the two methods:
library(SPCompute)

## Binary Covariate
parameter_binary_binary <- list(preva = 0.3, TraitSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_binary <- list(TraitMean = 0, TraitSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_binary_2 <- list(ResidualSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)


test_that("Check Power computations are consistent",{
  expect_equal(Compute_Power(parameter_binary_binary, response = "binary", covariate = "binary", n = 1000, method = "semi-sim"), Compute_Power(parameter_binary_binary, response = "binary", covariate = "binary", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_binary, response = "continuous", covariate = "binary", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_binary, response = "continuous", covariate = "binary", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_binary_2, response = "continuous", covariate = "binary", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_binary_2, response = "continuous", covariate = "binary", n = 1000, method = "expand"), tolerance = 0.01)
  }
)


## Continuous Covariate
parameter_binary_cont <- list(preva = 0.3, TraitSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_cont <- list(TraitMean = 0, TraitSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_cont_2 <- list(ResidualSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)


test_that("Check Power computations are consistent",{
  expect_equal(Compute_Power(parameter_binary_cont, response = "binary", covariate = "continuous", n = 1000, method = "semi-sim"), Compute_Power(parameter_binary_cont, response = "binary", covariate = "continuous", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_cont, response = "continuous", covariate = "continuous", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_cont, response = "continuous", covariate = "continuous", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_cont_2, response = "continuous", covariate = "continuous", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_cont_2, response = "continuous", covariate = "continuous", n = 1000, method = "expand"), tolerance = 0.01)
}
)

## No Covariate
parameter_binary_no <- list(preva = 0.3, TraitSD = 1, pG = 0.4, betaG = 1)
parameter_cont_no <- list(TraitMean = 0, TraitSD = 1, pG = 0.4, betaG = 1)
parameter_cont_no_2 <- list(ResidualSD = 1,  pG = 0.4, betaG = 1)

test_that("Check Power computations are consistent",{
  expect_equal(Compute_Power(parameter_binary_no, response = "binary", covariate = "none", n = 1000, method = "semi-sim"), Compute_Power(parameter_binary_no, response = "binary", covariate = "none", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_no, response = "continuous", covariate = "none", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_no, response = "continuous", covariate = "none", n = 1000, method = "expand"), tolerance = 0.01)
  expect_equal(Compute_Power(parameter_cont_no_2, response = "continuous", covariate = "none", n = 1000, method = "semi-sim"), Compute_Power(parameter_cont_no_2, response = "continuous", covariate = "none", n = 1000, method = "expand"), tolerance = 0.01)
}
)
