#' Compute the required power for binary response, a SNP G and two continuous covariates that are conditionally independent given G, using the Semi-Sim method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, muE, sigmaE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The power that can be achieved at the given sample size (using semi-sim method).
#' @noRd
Compute_Power_Sim_BCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, LargePowerApproxi = FALSE){
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (2*pG*qG + 2*pG^2)
    gamma02 <- muE2 - gammaG2 * (2*pG*qG + 2*pG^2)


    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}


    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)


    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)

      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)
    betaG <- parameters$betaG

    varG <- pG*qG
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_dom_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_dom_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    ### Simulate for SE: by averaging B times
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    set.seed(seed)

    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    pG <- parameters$pG^2
    qG <- 1 - pG

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)

    betaG <- parameters$betaG

    if((sigmaE1^2) <= (gammaG1^2) * qG*pG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * qG*pG)
    if((sigmaE2^2) <= (gammaG2^2) * qG*pG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * qG*pG)


    solveForbeta0_rec_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_rec_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    ### Simulate for SE: by averaging B times
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    set.seed(seed)

    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)

      X <- matrix(c(1,G,E1, E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
}



#' Compute the required power for binary response, a SNP G and two continuous covariates that are conditionally independent given G, using the empirical method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, muE, sigmaE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param seed An integer number that indicates the seed used for the simulation, by default is 123.
#' @return The power that can be achieved at the given sample size (computed from empirical power).
#' @noRd
Compute_Power_Emp_BCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8){
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (2*pG*qG + 2*pG^2)
    gamma02 <- muE2 - gammaG2 * (2*pG*qG + 2*pG^2)

    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "dominant"){
    preva <- parameters$preva
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)
    betaG <- parameters$betaG
    varG <- qG*pG
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_dom_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_dom_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B

  }
  else if(mode == "recessive") {
    preva <- parameters$preva
    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    pG <- parameters$pG^2
    qG <- 1 - pG

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)

    betaG <- parameters$betaG

    if((sigmaE1^2) <= (gammaG1^2) * qG*pG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * qG*pG)
    if((sigmaE2^2) <= (gammaG2^2) * qG*pG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * qG*pG)


    solveForbeta0_rec_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_rec_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  Power
}






#' Compute the required power for binary response, a SNP G and two binary covariates that are conditionally independent given G, using the Semi-Sim method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, pE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The power that can be achieved at the given sample size (using semi-sim method).
#' @noRd
Compute_Power_Sim_BBB <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  ComputeYgivenGE <- function(beta0,betaG, betaE1, betaE2, G, E1, E2, Y = 1){
    PYGE <- (exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2)^Y)/(1 + exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2))
    PYGE
  }
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)
      ComputeP <- function(beta0){
        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG^2)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG^2)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG^2)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG^2)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (2*qG*pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (2*qG*pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (2*qG*pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (2*qG*pG)

        P200 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 0) * (pG^2)
        P201 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 1) * (pG^2)
        P210 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 0) * (pG^2)
        P211 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 1) * (pG^2)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111 + P200 + P201 + P210 + P211
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)

      ComputeP <- function(beta0){

        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
    pG <- (parameters$pG)^2
    qG <- (1 - pG)
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)

      ComputeP <- function(beta0){

        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
}







#' Compute the required power for binary response, a SNP G and two binary covariates that are conditionally independent given G, using the empirical method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, muE, sigmaE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param seed An integer number that indicates the seed used for the simulation, by default is 123.
#' @return The power that can be achieved at the given sample size (computed from empirical power).
#' @noRd
Compute_Power_Emp_BBB <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8){
  correct <- c()
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  ComputeYgivenGE <- function(beta0,betaG, betaE1, betaE2, G, E1, E2, Y = 1){
    PYGE <- (exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2)^Y)/(1 + exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2))
    PYGE
  }
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)
      ComputeP <- function(beta0){
        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG^2)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG^2)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG^2)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG^2)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (2*qG*pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (2*qG*pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (2*qG*pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (2*qG*pG)

        P200 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 0) * (pG^2)
        P201 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 1) * (pG^2)
        P210 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 0) * (pG^2)
        P211 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 2, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 2, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 2, E = 1) * (pG^2)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111 + P200 + P201 + P210 + P211
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "dominant"){
    preva <- parameters$preva
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)

      ComputeP <- function(beta0){

        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "recessive") {
    preva <- parameters$preva
    pG <- (parameters$pG)^2
    qG <- (1 - pG)
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG
    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    solveForbeta0_add <- function(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2){
      qG <- 1 - pG
      gamma01 <- solveForgamma0(pE1,gammaG1, pG)
      gamma02 <- solveForgamma0(pE2,gammaG2, pG)

      ComputeP <- function(beta0){

        P000 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P001 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)
        P010 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 0) * (qG)
        P011 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 0, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 0, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 0, E = 1) * (qG)

        P100 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P101 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 0, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 0) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)
        P110 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 0) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 0) * (pG)
        P111 <- ComputeYgivenGE(beta0,betaG, betaE1, betaE2, G = 1, E1 = 1, E2 = 1) * ComputeEgivenG(gamma01,gammaG1,G = 1, E = 1) * ComputeEgivenG(gamma02,gammaG2,G = 1, E = 1) * (pG)

        P <- P000 + P001 + P010 + P011 + P100 + P101 + P110 + P111
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)
    beta0 <- solveForbeta0_add(preva, betaG, betaE1, betaE2, pG, pE1, pE2, gammaG1, gammaG2)
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  Power
}






#' Compute the required power for binary response, a SNP G and two covariate (one binary, one continuous) that are conditionally independent given G, using the Semi-Sim method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE and gammaG have to be vectors of length 2. The binary covariate is assumed to be the first covariate.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The power that can be achieved at the given sample size (using semi-sim method).
#' @noRd
Compute_Power_Sim_BBC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  ComputeYgivenGE <- function(beta0,betaG, betaE1, betaE2, G, E1, E2, Y = 1){
    PYGE <- (exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2)^Y)/(1 + exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2))
    PYGE
  }
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (2*pG*qG + 2*pG^2)


    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }

    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (pG)


    betaG <- parameters$betaG
    varG <- pG*qG

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG, pG))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
    pG <- parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }

    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (pG)


    betaG <- parameters$betaG
    varG <- pG*qG

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG, pG))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
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
}







#' Compute the required power for binary response, a SNP G and two covariate (one binary, one continuous) that are conditionally independent given G, using the empirical power.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE and gammaG have to be vectors of length 2. The binary covariate is assumed to be the first covariate.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param searchSizeGamma0 The interval radius for the numerical search of gamma0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The power that can be achieved at the given sample size (using semi-sim method).
#' @noRd
Compute_Power_Emp_BBC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  ComputeYgivenGE <- function(beta0,betaG, betaE1, betaE2, G, E1, E2, Y = 1){
    PYGE <- (exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2)^Y)/(1 + exp(beta0 + betaG * G + betaE1 * E1 + betaE2 * E2))
    PYGE
  }
  correct <- c()
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (2*pG*qG + 2*pG^2)

    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaE1 * E1 + betaE2 * E2 + betaG * G + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "dominant"){
    preva <- parameters$preva
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }

    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (pG)

    betaG <- parameters$betaG
    varG <- pG*qG

    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG, pG))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaE1 * E1 + betaE2 * E2 + betaG * G + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "recessive") {
    preva <- parameters$preva
    pG <- parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
          ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }

    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * (pG)

    betaG <- parameters$betaG
    varG <- pG*qG

    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG, pG))
        E1 <- gamma01 + gammaG1 * G + stats::rlogis(B)
        E1 <- ifelse(E1 >= 0, 1, 0)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaE1 * E1 + betaE2 * E2 + betaG * G + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
}
  Power
}



