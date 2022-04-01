## Testing for Power Computation using the two methods:
library(SPCompute)

## Binary Covariate
parameter_binary_binary <- list(preva = 0.3, TraitSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_binary <- list(TraitMean = 0, TraitSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_binary_2 <- list(ResidualSD = 1, gammaG = 0.1, pE = 0.4, pG = 0.4, betaG = 1, betaE = 1)


test_that("Parameters are checked correctly",{
  expect_equal(check_parameters(parameters = parameter_binary_binary, response = "binary", covariate = "binary"), TRUE)
  expect_equal(check_parameters(parameters = parameter_cont_binary, response = "continuous", covariate = "binary"), TRUE)
  expect_equal(check_parameters(parameters = parameter_cont_binary_2, response = "continuous", covariate = "binary"), TRUE)
}
)

test_that("Wrong parameters will be detected",{
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_binary, response = "binary", covariate = "continuous")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_binary, response = "continuous", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_binary, response = "continuous", covariate = "none")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_binary, response = "none", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_binary, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_binary_2, response = "binary", covariate = "binary")), FALSE)
}
)


## Continuous Covariate
parameter_binary_cont <- list(preva = 0.3, TraitSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_cont <- list(TraitMean = 0, TraitSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)
parameter_cont_cont_2 <- list(ResidualSD = 1, gammaG = 0.1, muE = 0.4, sigmaE = 1, pG = 0.4, betaG = 1, betaE = 1)



test_that("Parameters are checked correctly",{
  expect_equal(check_parameters(parameters = parameter_binary_cont, response = "binary", covariate = "continuous"), T)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont, response = "continuous", covariate = "continuous")), T)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont_2, response = "continuous", covariate = "continuous")), T)
}
)

test_that("Wrong parameters will be detected",{
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_cont, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont_2, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_cont, response = "continuous", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont, response = "continuous", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont_2, response = "continuous", covariate = "binary")), FALSE)
}
)


## No Covariate
parameter_binary_no <- list(preva = 0.3, TraitSD = 1, pG = 0.4, betaG = 1)
parameter_cont_no <- list(TraitMean = 0, TraitSD = 1, pG = 0.4, betaG = 1)
parameter_cont_no_2 <- list(ResidualSD = 1,  pG = 0.4, betaG = 1)

test_that("Parameters are checked correctly",{
  expect_equal(check_parameters(parameters = parameter_binary_cont, response = "binary", covariate = "none"), T)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont, response = "continuous", covariate = "none")), T)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_cont_2, response = "continuous", covariate = "none")), T)
}
)

test_that("Wrong parameters will be detected",{
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_no, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_no, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_no_2, response = "binary", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_binary_no, response = "continuous", covariate = "none")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_no, response = "continuous", covariate = "binary")), FALSE)
  expect_equal(suppressMessages(check_parameters(parameters = parameter_cont_no_2, response = "continuous", covariate = "binary")), FALSE)
}
)
