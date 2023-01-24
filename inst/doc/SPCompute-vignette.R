## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SPCompute)

## -----------------------------------------------------------------------------
para <- list(preva = 0.2, pG = 0.1, betaG = 0.1, betaE = 0.3, pE = 0.3, gammaG = 0)

## -----------------------------------------------------------------------------
Compute_Power(parameters = para, n = 2e4, covariate = "binary", mode = "additive", alpha = 0.05, method = "semi-sim")

## -----------------------------------------------------------------------------
Compute_Power(parameters = para, n = 2e4, covariate = "binary", mode = "additive", alpha = 0.05, method = "expand")

## -----------------------------------------------------------------------------
round(Compute_Size(parameters = para, PowerAim = 0.8, covariate = "binary", mode = "additive", alpha = 0.05, method = "semi-sim"))
round(Compute_Size(parameters = para, PowerAim = 0.8, covariate = "binary", mode = "additive", alpha = 0.05, method = "expand"))

## -----------------------------------------------------------------------------
round(Compute_Size(parameters = para, PowerAim = 0.8, covariate = "binary", mode = "dominant", alpha = 0.05, method = "semi-sim"))

round(Compute_Size(parameters = para, PowerAim = 0.8, covariate = "binary", mode = "dominant", alpha = 0.05, method = "expand"))

## -----------------------------------------------------------------------------
round(Compute_Size(parameters = para, PowerAim = 0.8, covariate = "binary", mode = "dominant", alpha = 0.05, method = "semi-sim", B = 5e5))

## -----------------------------------------------------------------------------
para <- list(TraitMean = 3, TraitSD = 1, pG = 0.1, betaG = 0.1, betaE = 0.3, pE = 0.3, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

## -----------------------------------------------------------------------------
para <- list(TraitMean = 3, ResidualSD = 1, pG = 0.1, betaG = 0.1, betaE = 0.3, pE = 0.3, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

## -----------------------------------------------------------------------------
para <- list(TraitMean = 3, ResidualSD = 1, pG = 0.1, betaG = 0.1, betaE = 0.3, pE = 0.3, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

para <- list(TraitMean = 30000, ResidualSD = 1, pG = 0.1, betaG = 0.1, betaE = 0.3, pE = 0.3, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

## -----------------------------------------------------------------------------
para <- list(TraitMean = 0, ResidualSD = 1, pG = 0.1, betaG = 0.1, betaE = 1000, pE = 0.5, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

para <- list(TraitMean = 0, ResidualSD = 1, pG = 0.1, betaG = 0.1, betaE = 0.001, pE = 0.1, gammaG = 0)
Compute_Power(parameters = para, n = 5e3, covariate = "binary", response = "continuous", mode = "additive", alpha = 0.05)
round(Compute_Size(parameters = para, PowerAim = 0.8, response = "continuous", mode = "additive", alpha = 0.05))

## -----------------------------------------------------------------------------
para <- list(preva = 0.2, pG = 0.1, betaG = 0.1, betaE = c(0.3, 0.2), 
             pE = 0.3, gammaG = c(0.15,0.25), gammaE = 0.1,
             muE = 0, sigmaE = sqrt(0.4))
Compute_Power_multi(parameters = para, n = 3000, mode = "additive",
                    covariate = c("binary", "continuous"), response = "binary")

