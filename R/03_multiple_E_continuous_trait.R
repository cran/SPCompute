#' Compute the required power for continuous response, a SNP G and two continuous covariates that are conditionally independent given G, using the Semi-Sim method.
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
Compute_Power_Sim_CCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, LargePowerApproxi = FALSE){
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
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
    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = SigmaErrorStage22)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B)*(1/(SigmaErrorStage1^2))
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
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG
    TraitMean <- parameters$TraitMean
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
    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = SigmaErrorStage22)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B)*(1/(SigmaErrorStage1^2))
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
    pG <- parameters$pG^2
    qG <- 1 - pG
    TraitMean <- parameters$TraitMean
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
    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = SigmaErrorStage22)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B)*(1/(SigmaErrorStage1^2))
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



#' Compute the required power for continuous response, a SNP G and two continuous covariates that are conditionally independent given G, using the empirical method.
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
Compute_Power_Emp_CCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8){
  correct <- c()
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
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
    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = SigmaErrorStage22)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n,sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y ~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "dominant"){
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG
    TraitMean <- parameters$TraitMean
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
    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = SigmaErrorStage22)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n,sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "recessive") {
    pG <- parameters$pG^2
    qG <- 1 - pG
    TraitMean <- parameters$TraitMean
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
    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)
    EGE1 <- gamma01 * EG + gammaG1 * EG2
    EGE2 <- gamma02 * EG + gammaG2 * EG2

    beta0 <- TraitMean - betaG * EG - betaE1 * muE1 - betaE2 * muE2
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    SigmaErrorStage21 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    SigmaErrorStage22 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)
    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      SigmaErrorStage1 <- sqrt((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)))
      if ((TraitSD^2 - (betaG^2) * varG - (betaE1^2) * (gammaG1^2)* varG - (betaE2^2) * (gammaG2^2)* varG - (betaE1^2)*(SigmaErrorStage21^2) - (betaE2^2)*(SigmaErrorStage22^2)) <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = SigmaErrorStage21)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = SigmaErrorStage22)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n,sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  Power
}







#' Compute the required power for continuous response, a SNP G and two binary covariates that are conditionally independent given G, using the Semi-Sim method.
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
Compute_Power_Sim_CBB <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
    pG <- parameters$pG
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
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
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (2*pG*qG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (2*pG*qG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG^2) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (2*pG*qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 2, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 2, E = 1)* (pG^2) - pE1*pE2


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2

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
      I <- I + X %*% t(X)
    }
    I <- I/B * (1/SigmaErrorStage1^2)
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
    TraitMean <- parameters$TraitMean
    qG <- (1 - parameters$pG)^2
    pG <- 1 - qG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (pG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (pG)

    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2

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
      I <- I + X %*% t(X)
    }
    I <- I/B * (1/SigmaErrorStage1^2)
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
    TraitMean <- parameters$TraitMean
    pG <- (parameters$pG)^2
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (pG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (pG)

    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2

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
      I <- I + X %*% t(X)
    }
    I <- I/B * (1/SigmaErrorStage1^2)
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




#' Compute the required power for continuous response, a SNP G and two binary covariates that are conditionally independent given G, using the empirical method.
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
Compute_Power_Emp_CBB <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  correct <- c()
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
    pG <- parameters$pG
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
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
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (2*pG*qG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (2*pG*qG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG^2) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (2*pG*qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 2, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 2, E = 1)* (pG^2) - pE1*pE2


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n = n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- mean(correct)
  }
  else if(mode == "dominant"){
    TraitMean <- parameters$TraitMean
    qG <- (1 - parameters$pG)^2
    pG <- 1 - qG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (pG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (pG)

    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2

    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n = n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- mean(correct)

  }
  else if(mode == "recessive") {
    TraitMean <- parameters$TraitMean
    pG <- (parameters$pG)^2
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    pE1 <- parameters$pE[1]
    qE1 <- 1 - pE1
    betaE1 <- parameters$betaE[1]
    gammaG2 <- parameters$gammaG[2]
    pE2 <- parameters$pE[2]
    qE2 <- 1 - pE2
    betaE2 <- parameters$betaE[2]
    betaG <- parameters$betaG

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    gamma01 <- solveForgamma0(pE1,gammaG1, pG)
    gamma02 <- solveForgamma0(pE2,gammaG2, pG)

    EG <- pG
    EG2 <- pG
    varG <- EG2 - (EG^2)

    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE1 * EG
    COV2G <- ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1) * (pG) - pE2 * EG
    COV12 <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 0, E = 1)* (qG) +
      ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * ComputeEgivenG(gamma0 = gamma02, gammaG = gammaG2, G = 1, E = 1)* (pG)

    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE1*qE1) - (betaE2^2)*(pE2*qE2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }


    beta0 <- TraitMean - betaG * EG - betaE1 * pE1 - betaE2 * pE2
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E2 <- gamma02 + gammaG2*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- ifelse(E2 >= 0, 1, 0)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n = n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    Power <- mean(correct)
  }
  Power
}










#' Compute the required power for continuous response, a SNP G and two covariate (one binary, one continuous) that are conditionally independent given G, using the Semi-Sim method.
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
Compute_Power_Sim_CBC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE

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

    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (2*pG*qG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG^2) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (2*pG*qG) +
      (gamma02 + gammaG2 * 2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 2, E = 1)) * (pG^2) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B) * (1/SigmaErrorStage1^2)
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
    TraitMean <- parameters$TraitMean
    qG <- (1 - parameters$pG)^2
    pG <- 1 - qG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    EG <- pG
    EG2 <- pG
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * EG
    betaG <- parameters$betaG
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (pG) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B) * (1/SigmaErrorStage1^2)
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
    TraitMean <- parameters$TraitMean
    pG <- (parameters$pG)^2
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE

    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    EG <- pG
    EG2 <- pG
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * EG
    betaG <- parameters$betaG
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (pG) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(1)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      I <- I + X %*% t(X)
    }
    I <- (I/B) * (1/SigmaErrorStage1^2)
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








#' Compute the required power for continuous response, a SNP G and two covariate (one binary, one continuous) that are conditionally independent given G, using the empirical power.
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
Compute_Power_Emp_CBC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, searchSizeGamma0 = 8, LargePowerApproxi = FALSE){
  ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
    PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
    PEG
  }
  correct <- c()
  if(mode == "additive"){
    TraitMean <- parameters$TraitMean
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE

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

    EG <- 2*pG*qG + 2* pG^2
    EG2 <- 2*pG*qG + 4* pG^2
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (2*pG*qG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG^2) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (2*pG*qG) +
      (gamma02 + gammaG2 * 2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 2, E = 1)) * (pG^2) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    mean(correct)
  }
  else if(mode == "dominant"){
    TraitMean <- parameters$TraitMean
    qG <- (1 - parameters$pG)^2
    pG <- 1 - qG

    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    EG <- pG
    EG2 <- pG
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * EG
    betaG <- parameters$betaG
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (pG) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    mean(correct)
  }
  else if(mode == "recessive"){
    TraitMean <- parameters$TraitMean
    pG <- (parameters$pG)^2
    qG <- 1 - pG
    gammaG1 <- parameters$gammaG[1]
    gammaG2 <- parameters$gammaG[2]
    betaE1 <- parameters$betaE[1]
    betaE2 <- parameters$betaE[2]

    pE <- parameters$pE
    qE <- 1 - pE

    muE <- parameters$muE
    sigmaE <- parameters$sigmaE

    solveForgamma0 <- function(pE,gammaG, pG){
      qG <- 1 - pG
      ComputePE <- function(gamma0){
        PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
        PE - pE
      }
      stats::uniroot(ComputePE, c(-searchSizeGamma0, searchSizeGamma0))$root
    }
    EG <- pG
    EG2 <- pG
    gamma01 <- solveForgamma0(pE, gammaG1, pG)
    gamma02 <- muE - gammaG2 * EG
    betaG <- parameters$betaG
    varG <- EG2 - (EG^2)

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    sigmaError <- sqrt(sigmaE^2 - (gammaG2^2) * varG)


    COV1G <- ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1) * (pG) - pE * EG
    COV2G <- (gamma02 + gammaG2*EG2) - muE * EG
    COV12 <- (gamma02) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 0, E = 1)) * (qG) +
      (gamma02 + gammaG2) * (ComputeEgivenG(gamma0 = gamma01, gammaG = gammaG1, G = 1, E = 1)) * (pG) - pE * muE


    if("TraitSD" %in% names(parameters)){
      TraitSD <- parameters$TraitSD
      ResidualVar <- TraitSD^2 - (betaE1^2)*(pE*qE) - (betaE2^2)*(sigmaE^2) - (betaG^2)*(varG) - 2*(betaG*betaE1*COV1G + betaG*betaE2*COV2G + betaE1*betaE2*COV12)
      SigmaErrorStage1 <- sqrt(ResidualVar)
      if (ResidualVar <= 0) {return(message("Error: TraitSD must be large enough to be compatible with other parameters"))}
    }
    else{
      SigmaErrorStage1 <- parameters$ResidualSD
    }

    beta0 <- TraitMean - betaE1 * pE - betaE2 * muE - betaG * EG
    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG,pG))
      E1 <- gamma01 + gammaG1*G + stats::rlogis(n)
      E1 <- ifelse(E1 >= 0, 1, 0)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rnorm(n, sd = SigmaErrorStage1)
      correct[i] <- summary(stats::glm(y~ G + E1 + E2, family = stats::gaussian()))$coefficients[2,4] <= alpha
    }
    mean(correct)
  }



}






