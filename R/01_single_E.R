
#' Check if the parameter list contains all the parameters required for the computation.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "preva" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "TraitSD" and "TraitMean" which represent the standard deviation and mean of the continuous trait.
#' #' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' #' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model.
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model.
#' @return TRUE or FALSE if all the parameters are correctly defined.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1),
#' betaE = log(1.1), muE = 0, sigmaE = 3, gammaG = log(2.1))
#'
#' SPCompute::check_parameters(parameters, "continuous", "continuous")
#' @export
check_parameters <- function(parameters, response, covariate){
  if(response == "binary"){
    if("continuous" %in% covariate){
      required_names <- c("preva", "pG", "betaG", "betaE", "muE", "sigmaE", "gammaG")
      if(any(! required_names %in% names(parameters))){
        message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
        return(F)
      }
      else{
        return(T)
      }
    }
    else if("binary" %in% covariate){
      required_names <- c("preva", "pG", "betaG", "betaE", "pE", "gammaG")
      if(any(! required_names %in% names(parameters))){
        message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
        return(F)
      }
      else{
        return(T)
      }
    }
    else{
      required_names <- c("preva", "pG", "betaG")
      if(any(! required_names %in% names(parameters))){
        message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
        return(F)
      }
      else{
        return(T)
      }
    }
  }
  if(response == "continuous"){
    if("continuous" %in% covariate){
      required_names <- c("TraitMean", "TraitSD", "pG", "betaG", "betaE", "muE", "sigmaE", "gammaG")
      required_names2 <- c("ResidualSD", "pG", "betaG", "betaE", "muE", "sigmaE", "gammaG")
      if(any(! required_names %in% names(parameters))){
        if(any(! required_names2 %in% names(parameters))){
          message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
          return(F)
        }
        else{
          T
        }
      }
      else{
        return(T)
      }
    }
    else if("binary" %in% covariate){
      required_names <- c("TraitMean", "TraitSD", "pG", "betaG", "betaE", "pE", "gammaG")
      required_names2 <- c("ResidualSD", "pG", "betaG", "betaE", "pE", "gammaG")
      if(any(! required_names %in% names(parameters))){
        if(any(! required_names2 %in% names(parameters))){
          message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
          return(F)
        }
        else{
          T
        }
      }
      else{
        return(T)
      }
    }
    else{
      required_names <- c("TraitMean", "TraitSD", "pG", "betaG")
      required_names2 <- c("ResidualSD", "pG", "betaG")
      if(any(! required_names %in% names(parameters))){
        if(any(! required_names2 %in% names(parameters))){
          message(paste("Error: There are missing parameters that need to be defined first:", paste(required_names[!required_names %in% names(parameters)], collapse = ",")))
          return(F)
        }
        else{
          T
        }
      }
      else{
        return(T)
      }
    }
  }
  else{
    message("Error: undefined response type")
    return(F)
  }
}






#' Compute the Power of an association study, at a given sample size, using the semi-simulation method
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @return The power that can be achieved at the given sample size.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1), betaE = log(1.1),
#' muE = 0, sigmaE = 3, gammaG = log(2.1))
#' SPCompute:::Compute_Power_Sim(parameters, n = 1000, B = 10000, "continuous", "continuous")
#' @noRd
Compute_Power_Sim <- function(parameters, n, B = 10000, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, searchSizeGamma0 = 8, searchSizeBeta0 = 8){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  if(response == "binary"){
    if(covariate == "binary"){
      ### Compute probability
      ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
        PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
        PEG
      }
      ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
        PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
        PYGE
      }
      ## Model 1:
      preva <- parameters$preva
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      ## Model 2:
      gammaG <- parameters$gammaG
      pE <- parameters$pE
      pG <- parameters$pG
      qG <- 1- pG
      if(mode == "additive"){
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_add <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG^2) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (2* pG * qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (2* pG * qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0, gammaG, G = 2)) * (pG^2)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_add(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "recessive"){
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
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_rec <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_rec(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B

        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }

      }
      else if(mode == "dominant"){
        qG <- qG^2
        pG <- 1- qG
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_dom <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_dom(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B

        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }

      }
      else{message("Error: mode must be either dominant, recessive or additive")}
    }
    else if(covariate == "continuous"){
      #### Computed parameters:
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (2*pG*qG + 2*pG^2)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

        solveForbeta0_add_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_con(preva, betaG, betaE, pG, gammaG)
        ### Simulate for SE: by averaging B times
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- pG*qG
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

        solveForbeta0_dom_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_con(preva, betaG, betaE, pG, gammaG)

        ### Simulate for SE: by averaging B times
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "recessive") {
        preva <- parameters$preva
        gammaG <- parameters$gammaG
        pG <- parameters$pG^2
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        if((sigmaE^2) <= (gammaG^2) * qG*pG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * qG*pG)

        solveForbeta0_rec_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_rec_con(preva, betaG, betaE, pG, gammaG)

        ### Simulate for SE: by averaging B times
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else if(covariate == "none"){
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_add_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <-  (qG^2)*exp(beta0)/(1 + exp(beta0)) + (2*pG*qG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (pG^2)*exp(beta0 + 2*betaG)/(1 + exp(beta0 + 2*betaG))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "recessive"){
        preva <- parameters$preva
        pG <- parameters$pG^2
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_rec_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }

        beta0 <- solveForbeta0_rec_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        qG <- (1 - parameters$pG)^2
        pG <- 1 - qG
        betaG <- parameters$betaG
        solveForbeta0_dom_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else{message("Error: covariate must be either continuous, binary or none")}
  }
  else if(response == "continuous"){
    if("TraitSD" %in% names(parameters)){
      if(covariate == "continuous"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaG^2) * VarG - (betaE^2) * (sigmaE^2) - 2*betaG*betaE*COVGE) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        if ((TraitSD^2 - (betaG^2) * VarG) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
        SigmaError1 <- sqrt(TraitSD^2 - (betaG^2) * VarG)
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
    else{
      if(covariate == "continuous"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
  }
  else(return(message("Error: undefined response type")))
}









#' Compute the sample size required for an association study, to achieve an aimed power, using the semi-simulation method.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param PowerAim An numeric value between 0 and 1 that indicates the aimed power of the study.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param upper.lim.n An integer number that indicates the largest sample size to be considered.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @return The sample size required to achieve the aimed power.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1),
#' betaE = log(1.1), muE = 0, sigmaE = 3, gammaG = log(2.1))
#' SPCompute:::Compute_Size_Sim(parameters, PowerAim = 0.8, B = 10000, "continuous", "continuous")
#' @noRd
Compute_Size_Sim <- function(parameters, PowerAim, B = 10000, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, upper.lim.n = 800000, searchSizeGamma0 = 8, searchSizeBeta0 = 8){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  if(response == "binary"){
    if(covariate == "binary"){
      ### Compute probability
      ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
        PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
        PEG
      }
      ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
        PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
        PYGE
      }
      ## Model 1:
      preva <- parameters$preva
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      ## Model 2:
      gammaG <- parameters$gammaG
      pE <- parameters$pE
      pG <- parameters$pG
      qG <- 1- pG
      if(mode == "additive"){
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_add <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG^2) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (2* pG * qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (2* pG * qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0, gammaG, G = 2)) * (pG^2)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_add(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) {return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else if(mode == "recessive"){
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
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_rec <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_rec(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B

        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }

      }
      else if(mode == "dominant"){
        qG <- qG^2
        pG <- 1- qG
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_dom <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_dom(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          E_lat <- gamma0 + gammaG*G + stats::rlogis(1)
          E <- ifelse(E_lat>=0, 1, 0)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B

        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }


      }
      else{message("Error: mode must be either dominant, recessive or additive")}
    }
    else if(covariate == "continuous"){
      #### Computed parameters:
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (2*pG*qG + 2*pG^2)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

        solveForbeta0_add_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_con(preva, betaG, betaE, pG, gammaG)
        ### Simulate for SE: by averaging B times
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          return(stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root)
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- pG*qG
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

        solveForbeta0_dom_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_con(preva, betaG, betaE, pG, gammaG)

        ### Simulate for SE: by averaging B times
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else if(mode == "recessive"){
        preva <- parameters$preva
        gammaG <- parameters$gammaG
        pG <- parameters$pG^2
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        if((sigmaE^2) <= (gammaG^2) * qG*pG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * qG*pG)

        solveForbeta0_rec_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(B, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(B)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_rec_con(preva, betaG, betaE, pG, gammaG)

        ### Simulate for SE: by averaging B times
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
          E <- gamma0 + gammaG*G + stats::rnorm(1,sd = sigmaError)
          X <- matrix(c(1,G,E), ncol = 1)
          eta <- beta0 + betaG*G + betaE*E
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else if(covariate == "none"){
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_add_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <-  (qG^2)*exp(beta0)/(1 + exp(beta0)) + (2*pG*qG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (pG^2)*exp(beta0 + 2*betaG)/(1 + exp(beta0 + 2*betaG))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)
        for (i in 1:B) {
          G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else if(mode == "recessive"){
        preva <- parameters$preva
        pG <- parameters$pG^2
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_rec_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }

        beta0 <- solveForbeta0_rec_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        qG <- (1 - parameters$pG)^2
        pG <- 1 - qG
        betaG <- parameters$betaG
        solveForbeta0_dom_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_non(preva, betaG, pG)
        I <- matrix(data = 0, nrow = 2, ncol = 2)
        set.seed(seed)

        for (i in 1:B) {
          G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
          X <- matrix(c(1,G), ncol = 1)
          eta <- beta0 + betaG*G
          weight <- stats::dlogis(eta)
          I <- I + weight* X %*% t(X)
        }
        I <- I/B
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else{message("Error: covariate must be either continuous, binary or none")}
  }
  else if(response == "continuous"){
    if("TraitSD" %in% names(parameters)){
      if(covariate == "continuous"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaG^2) * VarG - (betaE^2) * (sigmaE^2) - 2*betaG*betaE*COVGE) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        if ((TraitSD^2 - (betaG^2) * VarG) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
        SigmaError1 <- sqrt(TraitSD^2 - (betaG^2) * VarG)
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
    else{
      if(covariate == "continuous"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
          }
          else{
            compute_power_diff <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power - PowerAim
            }
            if(compute_power_diff(upper.lim.n) <=0) { return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))}
            stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          ((stats::qnorm(1-(alpha/2)) + stats::qnorm(PowerAim))^2)*(solve(I)[2,2])/(betaG^2)
        }
        else{
          compute_power_diff <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power - PowerAim
          }
          if(compute_power_diff(upper.lim.n) <=0) return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))
          stats::uniroot(compute_power_diff, c(0, upper.lim.n))$root
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
  }
  else(return(message("Error: undefined response type")))
}







#' Compute the Power of an association study, at a given sample size，using the idea of expanded dataset
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param n An integer number that indicates the sample size.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @return The power that can be achieved at the given sample size.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1), betaE = log(1.1),
#' muE = 0, sigmaE = 3, gammaG = log(2.1))
#' SPCompute:::Compute_Power_Expanded(parameters, n = 1000, "continuous", "continuous")
#' @noRd
Compute_Power_Expanded <- function(parameters, n, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, searchSizeGamma0 = 100, searchSizeBeta0 = 100){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  if(response == "binary"){
    if(covariate == "binary"){
      ### Compute probability
      ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
        PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
        PEG
      }
      ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
        PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
        PYGE
      }
      ## Model 1:
      preva <- parameters$preva
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      ## Model 2:
      gammaG <- parameters$gammaG
      pE <- parameters$pE
      pG <- parameters$pG
      qG <- 1- pG
      if(mode == "additive"){
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_add <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG^2) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (2* pG * qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (2* pG * qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0, gammaG, G = 2)) * (pG^2)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_add(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        count_vec <- c(floor(pG^2 * n), floor(2*qG*pG * n), (n - floor(pG^2 * n) - floor(2*qG*pG * n)), floor(floor(pG^2 * n) * ComputeEgivenG(gamma0,gammaG,2, E = 1)), (floor(pG^2 * n) - floor(floor(pG^2 * n) * ComputeEgivenG(gamma0,gammaG,2, E = 1))), floor(floor(2*qG*pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)), (floor(2*qG*pG * n) - floor(floor(2*qG*pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))), floor((n - floor(pG^2 * n) - floor(2*qG*pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1)), ((n - floor(pG^2 * n) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            floor(2*qG*pG * n)) - floor((n - floor(pG^2 * n) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           floor(2*qG*pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))))
        if(min(count_vec) <= 0){return("Error: sample size is too small to apply the data augmentation, try the alternative method.")}
        G <- c(rep(2,floor(pG^2 * n)), rep(1,floor(2*qG*pG * n)), rep(0, (n - floor(pG^2 * n) - floor(2*qG*pG * n))))
        E_2 <- c(rep(1,floor(floor(pG^2 * n) * ComputeEgivenG(gamma0,gammaG,2, E = 1))), rep(0, (floor(pG^2 * n) - floor(floor(pG^2 * n) * ComputeEgivenG(gamma0,gammaG,2, E = 1)))))
        E_1 <- c(rep(1,floor(floor(2*qG*pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))), rep(0, (floor(2*qG*pG * n) - floor(floor(2*qG*pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)))))
        E_0 <- c(rep(1,floor((n - floor(pG^2 * n) - floor(2*qG*pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))), rep(0, ((n - floor(pG^2 * n) - floor(2*qG*pG * n)) - floor((n - floor(pG^2 * n) - floor(2*qG*pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))) ))
        E <- c(E_2, E_1, E_0)
        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "recessive"){
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
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_rec <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_rec(preva, betaG, betaE, pG, pE, gammaG)
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        set.seed(seed)
        count_vec <- c(floor(pG * n), (n - floor(pG * n)), floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)), (floor(pG * n) - floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))),  floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))   ,  ((n - floor(pG * n)) - floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))))
        if(min(count_vec) <= 0){return("Error: sample size is too small to apply the data augmentation, try the alternative method.")}
        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))
        E_1 <- c(rep(1,floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))), rep(0, (floor(pG * n) - floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)))))
        E_0 <- c(rep(1,floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))), rep(0, ((n - floor(pG * n)) - floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))  )))
        E <- c(E_1, E_0)
        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }

      }
      else if(mode == "dominant"){
        qG <- qG^2
        pG <- 1- qG
        ### Solve for gamma0
        solveForgamma0 <- function(pE,gammaG, pG){
          qG <- 1 - pG
          ComputePE <- function(gamma0){
            PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
              ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
            PE - pE
          }
          stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
        }
        ### Solve for beta0
        solveForbeta0_dom <- function(preva, betaG, betaE, pG, pE, gammaG){
          qG <- 1 - pG
          gamma0 <- solveForgamma0(pE,gammaG, pG)
          ComputeP <- function(beta0){
            P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
              ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        beta0 <- solveForbeta0_dom(preva, betaG, betaE, pG, pE, gammaG)
        set.seed(seed)
        count_vec <- c(floor(pG * n), (n - floor(pG * n)), floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)), (floor(pG * n) - floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))),  floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))   ,  ((n - floor(pG * n)) - floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))))
        if(min(count_vec) <= 0){return("Error: sample size is too small to apply the data augmentation, try the alternative method.")}
        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))
        E_1 <- c(rep(1,floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1))), rep(0, (floor(pG * n) - floor(floor(pG * n) * ComputeEgivenG(gamma0,gammaG,1, E = 1)))))
        E_0 <- c(rep(1,floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))), rep(0, ((n - floor(pG * n)) - floor((n - floor(pG * n)) * ComputeEgivenG(gamma0,gammaG,0, E = 1))  )))
        E <- c(E_1, E_0)
        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }

      }
      else{message("Error: mode must be either dominant, recessive or additive")}
    }
    else if(covariate == "continuous"){
      ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
        PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
        PYGE
      }
      #### Computed parameters:
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (2*pG*qG + 2*pG^2)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)
        solveForbeta0_add_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1,2), size = 10000, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
            E <- gamma0 + gammaG * G + stats::rnorm(10000, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(10000)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_con(preva, betaG, betaE, pG, gammaG)

        if(floor(pG^2 * n) <= 1 | floor(2*qG*pG * n) <= 1 | (n - floor(pG^2 * n) - floor(2*qG*pG * n)) <= 1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}


        G <- c(rep(2,floor(pG^2 * n)), rep(1,floor(2*qG*pG * n)), rep(0, (n - floor(pG^2 * n) - floor(2*qG*pG * n))))
        muE2 <- gamma0 + gammaG * 2
        muE1 <- gamma0 + gammaG * 1
        muE0 <- gamma0

        E2 <- stats::qnorm((1:floor(pG^2 * n) - 0.375)/(floor(pG^2 * n) + 0.25), mean = muE2, sd = sigmaError)
        E1 <- stats::qnorm((1:floor(2*qG*pG * n) - 0.375)/(floor(2*qG*pG * n) + 0.25), mean = muE1, sd = sigmaError)
        E0 <- stats::qnorm((1:(n - floor(pG^2 * n) - floor(2*qG*pG * n)) - 0.375)/((n - floor(pG^2 * n) - floor(2*qG*pG * n)) + 0.25), mean = muE0, sd = sigmaError)

        E <- c(E2,E1,E0)

        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))

        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
        qG <- 1 - pG
        gammaG <- parameters$gammaG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        varG <- pG*qG
        I <- matrix(data = 0, nrow = 3, ncol = 3)
        if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

        solveForbeta0_dom_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = 10000, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(10000, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(10000)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_con(preva, betaG, betaE, pG, gammaG)

        if(floor(pG * n) <= 1 | (n - floor(pG * n)) <=1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}

        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))

        muE1 <- gamma0 + gammaG * 1
        muE0 <- gamma0


        E1 <- stats::qnorm((1:floor(pG * n) - 0.375)/(floor(pG * n) + 0.25), mean = muE1, sd = sigmaError)
        E0 <- stats::qnorm((1:(n - floor(pG * n)) - 0.375)/((n - floor(pG * n)) + 0.25), mean = muE0, sd = sigmaError)


        E <- c(E1,E0)

        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))

        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }

      }
      else if(mode == "recessive") {
        preva <- parameters$preva
        gammaG <- parameters$gammaG
        pG <- parameters$pG^2
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gamma0 <- muE - gammaG * (pG)
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        if((sigmaE^2) <= (gammaG^2) * qG*pG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
        sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * qG*pG)

        solveForbeta0_rec_con <- function(preva, betaG, betaE, pG, gammaG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            set.seed(seed)
            G <- sample(c(0,1), size = 10000, replace = TRUE, prob = c(qG,pG))
            E <- gamma0 + gammaG * G + stats::rnorm(10000, sd = sigmaError)
            y <- beta0 + betaG * G + betaE * E + stats::rlogis(10000)
            P <- mean(ifelse(y > 0, 1, 0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_rec_con(preva, betaG, betaE, pG, gammaG)

        if(floor(pG * n) <= 1 | (n - floor(pG * n)) <=1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}


        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))

        muE1 <- gamma0 + gammaG * 1
        muE0 <- gamma0

        E1 <- stats::qnorm((1:floor(pG * n) - 0.375)/(floor(pG * n) + 0.25), mean = muE1, sd = sigmaError)
        E0 <- stats::qnorm((1:(n - floor(pG * n)) - 0.375)/((n - floor(pG * n)) + 0.25), mean = muE0, sd = sigmaError)

        E <- c(E1,E0)

        data <- data.frame(G = c(G,G), E = c(E,E))
        ## Attach Y into the dataframe, and compute the weight of each observation
        data$Y <- c(rep(1,n), rep(0,n))
        data$weights <- c(ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 1), ComputeYgivenGE(beta0, betaG, betaE, G = G, E = E, Y = 0))
        mod <- suppressWarnings(stats::glm(Y~G+E, family = stats::binomial(link = 'logit'), weights = data$weights, data = data))

        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else if(covariate == "none"){
      ComputeYgivenG <- function(beta0,betaG, G, Y = 1){
        PYG <- (exp(beta0 + betaG * G)^Y)/(1 + exp(beta0 + betaG * G))
        PYG
      }
      if(mode == "additive"){
        preva <- parameters$preva
        pG <- parameters$pG
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_add_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <-  (qG^2)*exp(beta0)/(1 + exp(beta0)) + (2*pG*qG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (pG^2)*exp(beta0 + 2*betaG)/(1 + exp(beta0 + 2*betaG))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_add_non(preva, betaG, pG)
        if(floor(pG^2 * n) <= 1 | floor(2*qG*pG * n) <= 1 | (n - floor(pG^2 * n) - floor(2*qG*pG * n)) <= 1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}

        G <- c(rep(2,floor(pG^2 * n)), rep(1,floor(2*qG*pG * n)), rep(0, (n - floor(pG^2 * n) - floor(2*qG*pG * n))))
        data <- data.frame(G = c(G,G), Y = c(rep(1, n), rep(0, n)), weights = c(ComputeYgivenG(beta0, betaG, G, Y = 1),ComputeYgivenG(beta0, betaG, G, Y = 0)))
        mod <- suppressWarnings(stats::glm(Y~G, family = stats::binomial(), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE <- summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "recessive"){
        preva <- parameters$preva
        pG <- parameters$pG^2
        qG <- 1 - pG
        betaG <- parameters$betaG
        solveForbeta0_rec_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }

        beta0 <- solveForbeta0_rec_non(preva, betaG, pG)
        if(floor(pG * n) <= 1 | (n - floor(pG * n)) <=1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}

        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))
        data <- data.frame(G = c(G,G), Y = c(rep(1, n), rep(0, n)), weights = c(ComputeYgivenG(beta0, betaG, G, Y = 1),ComputeYgivenG(beta0, betaG, G, Y = 0)))
        mod <- suppressWarnings(stats::glm(Y~G, family = stats::binomial(), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE <- summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else if(mode == "dominant"){
        preva <- parameters$preva
        qG <- (1 - parameters$pG)^2
        pG <- 1 - qG
        betaG <- parameters$betaG
        solveForbeta0_dom_non <- function(preva, betaG, pG){
          qG <- 1 - pG
          ComputeP <- function(beta0){
            P <- (pG)*exp(beta0 + betaG)/(1 + exp(beta0 + betaG)) + (qG)*exp(beta0)/(1 + exp(beta0))
            P - preva
          }
          stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
        }
        beta0 <- solveForbeta0_dom_non(preva, betaG, pG)
        if(floor(pG * n) <= 1 | (n - floor(pG * n)) <=1) {return(message("The sample size is not large enough to apply this method. Consider the alternative method based on semi-simulation."))}

        G <- c(rep(1,floor(pG * n)), rep(0, (n - floor(pG * n))))
        data <- data.frame(G = c(G,G), Y = c(rep(1, n), rep(0, n)), weights = c(ComputeYgivenG(beta0, betaG, G, Y = 1),ComputeYgivenG(beta0, betaG, G, Y = 0)))
        mod <- suppressWarnings(stats::glm(Y~G, family = stats::binomial(), weights = data$weights, data = data))
        if(LargePowerApproxi){
          SE <- summary(mod)$coefficients[2,2]
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE <- summary(mod)$coefficients[2,2]
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{return(message("Error: mode must be either dominant, recessive or additive"))}
    }
    else{message("Error: covariate must be either continuous, binary or none")}
  }
  else if(response == "continuous"){
    if("TraitSD" %in% names(parameters)){
      if(covariate == "continuous"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaG^2) * VarG - (betaE^2) * (sigmaE^2) - 2*betaG*betaE*COVGE) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))
          if (sigmaE^2 - (gammaG^2) * VarG <= 0) {return(message("Error: sigmaE must be large enough to be compatible with other parameters"))}
          SigmaError2 <- sqrt(sigmaE^2 - (gammaG^2) * VarG)
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          if ((TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          beta0 <- TraitMean - betaG * EG - betaE * muE

          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          if ((TraitSD^2 - (betaE^2)*(sigmaE^2) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
          SigmaError1 <- sqrt(TraitSD^2 - (betaE^2)*(pE*qE) - (betaG^2) * VarG - (2 * betaG * betaE * COVGE))

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        TraitMean <- parameters$TraitMean
        TraitSD <- parameters$TraitSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        if ((TraitSD^2 - (betaG^2) * VarG) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
        SigmaError1 <- sqrt(TraitSD^2 - (betaG^2) * VarG)
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
    else{
      if(covariate == "continuous"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        muE <- parameters$muE
        sigmaE <- parameters$sigmaE
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          gamma0 <- muE - gammaG*EG
          EGE <- gamma0 * EG + gammaG * EG2
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "binary"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        betaE <- parameters$betaE
        pG <- parameters$pG
        qG <- 1 - pG
        pE <- parameters$pE
        qE <- 1 - pE
        muE <- pE
        sigmaE <- sqrt(pE * qE)
        gammaG <- parameters$gammaG

        if(mode == "additive"){
          EG <- 2*pG*qG + 2* pG^2
          EG2 <- 2*pG*qG + 4* pG^2
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG^2) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (2*pG*qG) +
              exp(gamma0 + gammaG * 2)/(1+ exp(gamma0 + gammaG * 2)) * (pG^2)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (2*pG*qG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1)) + 2 * (pG^2) * exp(gamma0 + gammaG * 2)/(1+exp(gamma0 + gammaG * 2))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "dominant"){
          qG <- qG^2
          pG <- 1 - qG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE

          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <-  Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else if(mode == "recessive"){
          pG <- pG^2
          qG <- 1 - pG
          EG <- pG
          EG2 <- pG
          VarG <- EG2 - EG^2
          FunctionForGamma0 <- function(gamma0){
            EE <- exp(gamma0 + gammaG * 0)/(1+ exp(gamma0 + gammaG * 0)) * (qG) + exp(gamma0 + gammaG * 1)/(1+ exp(gamma0 + gammaG * 1)) * (pG)
            EE - muE
          }
          gamma0 <- stats::uniroot(FunctionForGamma0, c(-8, 8))$root
          EGE <- (pG)* 1 * exp(gamma0 + gammaG * 1)/(1+exp(gamma0 + gammaG * 1))
          COVGE <- EGE - EG*muE
          EE2 <- muE^2 + sigmaE^2
          I <- matrix(1, nrow = 3, ncol = 3)
          I[1,2] <- EG
          I[1,3] <- muE
          I[2,2] <- EG2
          I[2,3] <- EGE
          I[3,3] <- EE2
          I <- Matrix::forceSymmetric(I/(SigmaError1^2), uplo = "U")
          if(LargePowerApproxi){
            SE <- sqrt((solve(I)[2,2]))/sqrt(n)
            return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
          }
          else{
            compute_power <- function(n){
              ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
              SE = sqrt((solve(I)[2,2]))/sqrt(n)
              ### Once know this SE of betaG hat, compute its power at this given sample size n:
              Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
              Power
            }
            compute_power(n)
          }
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
      }
      else if(covariate == "none"){
        SigmaError1 <- parameters$ResidualSD
        betaG <- parameters$betaG
        pG <- parameters$pG
        qG <- 1 - pG
        if(mode == "additive"){
          I <- matrix(1, nrow = 2, ncol = 2)
          I[1,2] <- 2*pG*qG + 2*pG^2
          I[2,1] <- 2*pG*qG + 2*pG^2
          I[2,2] <- 2*pG*qG + 4*pG^2
          VarG <- 2*pG*qG + 4*pG^2 - (2*pG*qG + 2*(pG^2))^2
        }
        else if(mode == "dominant"){
          p <- 1 - qG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else if(mode == "recessive"){
          p <- pG^2
          I <- matrix(p, nrow = 2, ncol = 2)
          I[1,1] <- 1
          VarG <- p * (1-p)
        }
        else(return(message("Error: mode must be either dominant, recessive or additive")))
        I <- I/(SigmaError1^2)
        if(LargePowerApproxi){
          SE <- sqrt((solve(I)[2,2]))/sqrt(n)
          return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
        }
        else{
          compute_power <- function(n){
            ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
            SE = sqrt((solve(I)[2,2]))/sqrt(n)
            ### Once know this SE of betaG hat, compute its power at this given sample size n:
            Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
            Power
          }
          compute_power(n)
        }
      }
      else{message("Error: covariate must be either continuous, binary or none")}
    }
  }
  else(return(message("Error: undefined response type")))
}






#' Compute the sample size of a given study using the method of representative dataset
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param PowerAim An numeric value between 0 and 1 that indicates the aimed power of the study.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param lower.lim.n An integer number that indicates the smallest sample size to be considered, only for "expand" method.
#' @param upper.lim.n An integer number that indicates the largest sample size to be considered.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The sample size required to achieve the aimed power.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1), betaE = log(1.1),
#' muE = 0, sigmaE = 3, gammaG = log(2.1))
#' SPCompute:::Compute_Size_Expanded(parameters, PowerAim = 0.8, response = "continuous",
#' covariate = "continuous", mode = "additive")
#' @noRd
Compute_Size_Expanded <- function(parameters, PowerAim, response, covariate, mode, alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, lower.lim.n = 1000, upper.lim.n = 800000, searchSizeGamma0 = 100, searchSizeBeta0 = 100){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  compute_power_diff <- function(n){
    ### Once know this SE of betaG hat, compute its power at this given sample size n:
    Power = Compute_Power_Expanded(parameters = parameters, n = n, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi)
    Power - PowerAim
  }
  if(is.numeric(Compute_Power_Expanded(parameters = parameters, n = lower.lim.n, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi))){
    if(compute_power_diff(lower.lim.n) >= 0) return(message(paste("The required sample size is estimated to be smaller than lower.lim.n, which is", lower.lim.n, sep = " ")))
    if(compute_power_diff(upper.lim.n) <=0) return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))
    stats::uniroot(compute_power_diff, c(lower.lim.n, upper.lim.n))$root
  }
  else{
    message("Error with data augmentation occurs: consider trying larger lower.lim.n or other search limits")
    return(Compute_Power_Expanded(parameters = parameters, n = lower.lim.n, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi))
  }
}




#' Compute the sample size required for an association study, to achieve an aimed power, using the semi-simulation method (Implementing in a simpler but slower way).
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param PowerAim An numeric value between 0 and 1 that indicates the aimed power of the study.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param upper.lim.n An integer number that indicates the largest sample size to be considered.
#' @return The sample size required to achieve the aimed power.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1),
#' betaE = log(1.1), muE = 0, sigmaE = 3, gammaG = log(2.1))
#' SPCompute:::Compute_Size_Sim_simpler(parameters, PowerAim = 0.8,
#' B = 10000, "continuous", "continuous")
#' @noRd
Compute_Size_Sim_simpler <- function(parameters, PowerAim, B = 10000, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, upper.lim.n = 800000){
  if(check_parameters(parameters, response, covariate) != TRUE){return(message("Define the above missing parameters before continuing"))}
  compute_power_diff <- function(n){
    ### Once know this SE of betaG hat, compute its power at this given sample size n:
    Power = Compute_Power_Sim(parameters, n, B = B, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi)
    Power - PowerAim
  }
  if(compute_power_diff(upper.lim.n) <=0) return(message(paste("The required sample size is estimated to be larger than upper.lim.n, which is", upper.lim.n, sep = " ")))
  stats::uniroot(compute_power_diff, c(1, upper.lim.n))$root
}









#' Compute the Power of an association study, at a given sample size.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param method An character that is either "semi-sim" (default) or "expand" using the idea of representative dataset. This specifies the method being used to compute the power/sample size when the trait is binary using logistic regression. The default method will be faster for large sample size computation.
#' @param n An integer number that indicates the sample size.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
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
#' betaE = log(1.1), muE = 0, sigmaE = 3, gammaG = log(2.1))
#'
#' Compute_Power(parameters, n = 1000, response = "continuous",
#' covariate = "continuous", method = "semi-sim")

#' @export
Compute_Power <- function(parameters, n, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, searchSizeGamma0 = 100, searchSizeBeta0 = 100, B = 10000, method = "semi-sim"){
  if(method == "semi-sim"){
    Compute_Power_Sim(parameters, n, B = B, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi, searchSizeGamma0 = searchSizeGamma0, searchSizeBeta0 = searchSizeBeta0)
  }
  else if (method == "expand"){
    Compute_Power_Expanded(parameters, n, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi, searchSizeGamma0 = searchSizeGamma0, searchSizeBeta0 = searchSizeBeta0)
  }
  else{
    return(message("Error: The selected method is not defined."))
  }
}






#' Compute the sample size of an association study, to achieve a target power.
#'
#' @param parameters A list of parameters that contains all the required parameters in the model. If response is "binary", this list needs to contain "prev" which denotes the prevalence of the disease (or case to control ratio for case-control sampling). If response is continuous, the list needs to contain "traitSD" and "traitMean" which represent the standard deviation and mean of the continuous trait.
#' If covariate is not "none", a parameter "gammaG" needs to be defined to capture the dependence between the SNP and the covariate (through linear regression model if covariate is continuous, and logistic model if covariate is binary). If covariate is "binary", list needs to contains "pE" that defines the frequency of the covariate. If it is continuous, list needs to contain "muE" and "sigmaE" to define
#' its mean and standard deviation. The MAF is defined as "pG", with HWE assumed to hold.
#' @param method An character that is either "semi-sim" (default) or "expand" using the idea of representative dataset. This specifies the method being used to compute the power/sample size when the trait is binary using logistic regression. The default method will be faster for large sample size computation.
#' @param PowerAim An numeric value between 0 and 1 that indicates the aimed power of the study.
#' @param response A string of either "binary" or "continuous", indicating the type of response/trait variable in the model, by default is "binary"
#' @param covariate A string of either "binary", "continuous" or "none" indicating the type of covariate E in the model, by default is "binary".
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000 (Should only be changed if computation uses semi-simulation method).
#' @param lower.lim.n An integer number that indicates the smallest sample size to be considered, only for "expand" method.
#' @param upper.lim.n An integer number that indicates the largest sample size to be considered.
#' @return The required sample size.
#' @examples
#' parameters <- list(TraitMean = 0.3, TraitSD = 1, pG = 0.2, betaG = log(1.1),
#' betaE = log(1.1), muE = 0, sigmaE = 3, gammaG = log(2.1))
#'
#' Compute_Size(parameters, PowerAim = 0.8, response = "continuous",
#' covariate = "continuous", method = "semi-sim")
#' @export
Compute_Size <- function(parameters, PowerAim, response = "binary", covariate = "binary", mode = "additive", alpha = 0.05, seed = 123, LargePowerApproxi = FALSE, searchSizeGamma0 = 100, searchSizeBeta0 = 100, B = 10000, method = "semi-sim", lower.lim.n = 1000, upper.lim.n = 800000){
  if(method == "semi-sim"){
    Compute_Size_Sim(parameters = parameters, PowerAim = PowerAim, B = B, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi, searchSizeGamma0 = searchSizeGamma0, searchSizeBeta0 = searchSizeBeta0, upper.lim.n = upper.lim.n)
  }
  else if (method == "expand"){
    Compute_Size_Expanded(parameters = parameters, PowerAim = PowerAim, response = response, covariate = covariate, mode = mode, alpha = alpha, seed = seed, LargePowerApproxi = LargePowerApproxi, searchSizeGamma0 = searchSizeGamma0, searchSizeBeta0 = searchSizeBeta0, upper.lim.n = upper.lim.n, lower.lim.n = lower.lim.n)
  }
  else{
    return(message("Error: The selected method is not defined."))
  }
}













