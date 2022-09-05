#' Convert the prevalence value to the intercept value beta0.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param seed An integer number that indicates the seed used for the simulation if needed, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param B An integer number that indicates the number of simulated sample to use if needed, by default is 10000.
#' @return The corresponding gamma0, beta0 and residual variance of E (if applicable).
#' @examples
#' convert_preva_to_intercept(parameters = list(preva = 0.2, betaG = 0.6, betaE = c(0.9),
#' gammaG = c(0.2), muE = c(0), sigmaE = c(1), pG = 0.3), covariate = "continuous")
#' @export
convert_preva_to_intercept <- function(parameters, mode = "additive", covariate = "binary", seed = 123, B = 10000, searchSizeBeta0 = 8, searchSizeGamma0 = 8){
  check_parameters(parameters = parameters, response = "binary", covariate = covariate)
  preva <- parameters$preva
  betaG <- parameters$betaG
  if(covariate == "binary"){
    ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
      PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
      PEG
    }
    ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
      PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
      PYGE
    }

    pG <- parameters$pG
    pE <- parameters$pE
    betaE <- parameters$betaE
    gammaG <- parameters$gammaG

    if(mode == "recessive"){
      pG <- pG^2
      qG <- 1- pG
      ### Solve for gamma0
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
          PE - pE
        }
        stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0_11 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0_11(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))
    }
    else if(mode == "dominant"){
      qG <- (1-pG)^2
      pG <- 1 - qG
      ### Solve for gamma0
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
          PE - pE
        }
        stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0_12 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0_12(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))
    }
    else{
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
          PE - pE
        }
        stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0_13 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG^2) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (2* pG * qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (2* pG * qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0, gammaG, G = 2)) * (pG^2)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0_13(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))

    }
  }
  else if(covariate == "continuous"){
    pG <- parameters$pG
    qG <- 1 - pG
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE
    betaE <- parameters$betaE
    gammaG <- parameters$gammaG

    if(mode == "additive"){
      gamma0 <- muE - gammaG * (2*pG*qG + 2*pG^2)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)
      solveForbeta0_21 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1,2), size = B, replace = T, prob = c(qG^2, 2*pG*qG, pG^2))
          E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0_21(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }
    else if(mode == "dominant"){
      pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
      qG <- 1 - pG
      gammaG <- parameters$gammaG
      muE <- parameters$muE
      sigmaE <- parameters$sigmaE
      gamma0 <- muE - gammaG * (pG)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- pG*qG
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)
      solveForbeta0_22 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1), size = B, replace = T, prob = c(qG,pG))
          E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0_22(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }
    else{
      pG <- parameters$pG^2
      qG <- 1 - pG
      gammaG <- parameters$gammaG
      muE <- parameters$muE
      sigmaE <- parameters$sigmaE
      gamma0 <- muE - gammaG * (pG)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- pG*qG
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

      solveForbeta0_23 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1), size = B, replace = T, prob = c(qG,pG))
          E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0_23(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }

  }
  else{
    ComputeYgivenG <- function(beta0,betaG, G, Y = 1){
      PYG <- (exp(beta0 + betaG * G)^Y)/(1 + exp(beta0 + betaG * G))
      PYG
    }
    if(mode == "dominant"){
      solveForbeta0_32 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (pG)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      pG <- (parameters$pG^2) + 2*((1-parameters$pG) * parameters$pG)
      beta0 <- solveForbeta0_32(preva, betaG, pG)
      return(beta0)
    }
    else if(mode == "additive"){
      solveForbeta0_31 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG^2) + ComputeYgivenG(beta0,betaG, G = 2, Y = 1) * (pG^2) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (2*pG*qG)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0_31(preva, betaG, parameters$pG)
      return(beta0)
    }
    else{
      solveForbeta0_33 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (pG)
          P - preva
        }
        stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0), tol = (.Machine$double.eps^0.5) )$root
      }
      pG <- (parameters$pG^2)
      beta0 <- solveForbeta0_33(preva, betaG, pG)
      return(beta0)
    }

  }
}











#' Compute the Power of an association study at a given sample size, accommodating more than one covariates, using the Semi-Simulation method.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param n An integer number that indicates the sample size.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A vector of length two with elements being either c("binary", "continuous"),c("binary", "binary") or c("continuous", "continuous"),  indicating the type of covariate E in the model.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000 (Should only be changed if computation uses semi-simulation method).
#' @return The power that can be achieved at the given sample size.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1),
#' betaE = c(log(1.1), log(1.2)),
#' muE = 0, sigmaE = 3, gammaG = c(log(2.1), log(2.2)), pE = 0.4)
#' SPCompute::Compute_Power_multi(parameters, n = 1000, response = "continuous",
#' covariate = c("binary", "continuous"))
#' @export
Compute_Power_multi <- function(parameters, n, response = "binary", covariate, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE, B = 10000){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  if(response == "binary"){
    if(all(covariate == c("binary", "binary"))){
      Compute_Power_Sim_BBB(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, searchSizeGamma0 = searchSizeGamma0, LargePowerApproxi = LargePowerApproxi)
    }
    else if(all(covariate == c("binary", "continuous"))){
      Compute_Power_Sim_BBC(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, searchSizeGamma0 = searchSizeGamma0, LargePowerApproxi = LargePowerApproxi)
    }
    else if(all(covariate == c("continuous", "continuous"))){
      Compute_Power_Sim_BCC(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, LargePowerApproxi = LargePowerApproxi)
    }
    else{
      return(message("covariate only supports c(\"binary\", \"binary\"), c(\"binary\", or \"continuous\"), c(\"continuous\", \"continuous\")"))
    }
  }
  else if(response == "continuous") {
    if(all(covariate == c("binary", "binary"))){
      Compute_Power_Sim_CBB(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, searchSizeGamma0 = searchSizeGamma0, LargePowerApproxi = LargePowerApproxi)
    }
    else if(all(covariate == c("binary", "continuous"))){
      Compute_Power_Sim_CBC(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, searchSizeGamma0 = searchSizeGamma0, LargePowerApproxi = LargePowerApproxi)
    }
    else if(all(covariate == c("continuous", "continuous"))){
      Compute_Power_Sim_CCC(n = n, B = B, parameters = parameters, mode = mode, alpha = alpha, seed = seed, searchSizeBeta0 = searchSizeBeta0, LargePowerApproxi = LargePowerApproxi)
    }
    else{
      return(message("covariate only supports c(\"binary\", \"binary\"), c(\"binary\", or \"continuous\"), c(\"continuous\", \"continuous\")"))
    }
  }
  else{
    return(message("Response only allows 'binary' or 'continuous'."))
  }
}









#' Compute the sample size of an association study to achieve a target power for multiple E's, using semi-sim.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param PowerAim An numeric value between 0 and 1 that indicates the aimed power of the study.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate Same as in Compute_Power_multi.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param upper.lim.n An integer number that indicates the largest sample size to be considered.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000 (Should only be changed if computation uses semi-simulation method).
#' @return The required sample size.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2,
#' betaG = log(1.1), betaE = c(log(1.1), log(1.2)),
#' muE = 0, sigmaE = 3, gammaG = c(log(2.1), log(2.2)), pE = 0.4)
#' SPCompute::Compute_Size_multi(parameters, PowerAim = 0.8,
#' response = "continuous", covariate = c("binary", "continuous"))
#' @export
Compute_Size_multi <- function(parameters, PowerAim, response = "binary", covariate, mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, searchSizeGamma0 = 100, searchSizeBeta0 = 100, B = 10000, upper.lim.n = 800000){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  compute_power_diff <- function(n){
    ### Once know this SE of betaG hat, compute its power at this given sample size n:
    Power = Compute_Power_multi(parameters = parameters, n = n, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi, searchSizeGamma0 = searchSizeGamma0, searchSizeBeta0 = searchSizeBeta0, B = B)
    Power - PowerAim
  }
  if(compute_power_diff(upper.lim.n) <=0) return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))
  round(stats::uniroot(compute_power_diff, c(1, upper.lim.n))$root,0)
}








