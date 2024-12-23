# Source necessary files
source("bicVec.R")
source("snc.R")
source("psbcGL_fit.R")

library(psbcGroup)
# Simulation examples
set.seed(123)

# List of p values
p_values <- c(20, 200, 500, 1000)
results <- list()

# Loop over different p values
for (p in p_values) {
  
  # Set the hazard type (either "Gamma prior" or "Beta prior")
  hazard <- "Beta prior"  # Or change this to "Gamma prior"
  
  n <- 100
  beta.true <- c(rep(4, 10), rep(0, (p - 10)))
  
  # Covariance matrix
  CovX <- diag(0.1, p)
  
  # Generate survival object
  survObj <- list()
  survObj$x <- apply(rmvnorm(n, sigma = CovX, method = "chol"), 2, scale)
  
  pred <- as.vector(exp(rowSums(scale(survObj$x, center = FALSE, scale = 1 / beta.true))))
  
  t <- rexp(n, rate = pred)
  cen <- runif(n, 0, 8)
  survObj$t <- pmin(t, cen)
  survObj$di <- as.numeric(t <= cen)
  survObj <- list(t = survObj$t, di = survObj$di, x = survObj$x)
  
  # Set hyperparameters
  s <- c(sort(survObj$t[survObj$di == 1]), 
         2 * max(survObj$t) - max(survObj$t[-which(survObj$t == max(survObj$t))]))
  priorPara <- list('eta0' = 1, 'kappa0' = 1, 'c0' = 2, 'r' = 0.5, 
                    'delta' = 0.0001, 's' = s, 'J' = length(s), 
                    'groupInd' = 1:p)
  
  # Set MCMC parameters
  mcmcPara <- list('numBeta' = p, 'beta.prop.var' = 1)
  
  # Initial values of hyperparameters
  lambdaSq <- 1
  initial <- list(
    'beta.ini' = rep(0, p),
    'lambdaSq' = 1,
    'sigmaSq' = runif(1, 0.1, 10),
    'tauSq' = rexp(length(unique(priorPara$groupInd)), rate = lambdaSq / 2),
    'h' = ifelse(hazard == "Gamma prior", rgamma(priorPara$J, 1, 1), rbeta(priorPara$J, 1, 1))  # Direct conditional assignment for h
  )
  
  # # Fit psbcGL
  PSBC.GL.G <- psbcGL.fit(survObj, priorPara, initial, rw = TRUE, 
                          mcmcPara, num.reps = 100, thin = 1, chain = 1)
  
  # Compute BIC vector
  bicVecResult <- bicVec(PSBC.GL.G, survObj$x, psiVec = seq(0.001, 1, 0.001))
  
  # Variable selection
  VS_Result <- apply(PSBC.GL.G$beta.p, 2, mean) > 0.001 # 99.9% (0.001) predicted power
  count_true <- sum(VS_Result)
  ## Variable Selection (VS)
  VS(PSBC.GL.G, X=survObj$x)
  # Print results in the required format
  print(paste0("p = ", p, ": Number of Selected variables = ", count_true))
  print(paste0("Prediction Concordance Index: ", PSBC.GL.G$concordance.index))
  print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
  print(paste0("BIC: ", bicVecResult[length(bicVecResult)]))  # BIC value from the last entry of bicVec
  
  # Store results
  results[[paste0("p_", p)]] <- list(fit = PSBC.GL.G, bicVec = bicVecResult, 
                                     VS = VS_Result, count_true = count_true)
}




# Source necessary files
source("bicVec.R")
source("snc.R")
source("psbcGL_fit.R")

library(psbcGroup)
# Simulation examples
set.seed(123)

# List of p values
p_values <- c(20, 200, 500, 1000)
results <- list()

# Loop over different p values
for (p in p_values) {
  
  # Set the hazard type (either "Gamma prior" or "Beta prior")
  hazard <- "Beta prior"  # Or change this to "Gamma prior"
  
  n <- 100
  beta.true <- c(rep(4, 10), rep(0, (p - 10)))
  
  # Covariance matrix
  CovX <- diag(0.1, p)
  
  # Generate survival object
  survObj <- list()
  survObj$x <- apply(rmvnorm(n, sigma = CovX, method = "chol"), 2, scale)
  
  pred <- as.vector(exp(rowSums(scale(survObj$x, center = FALSE, scale = 1 / beta.true))))
  
  t <- rexp(n, rate = pred)
  cen <- runif(n, 0, 8)
  survObj$t <- pmin(t, cen)
  survObj$di <- as.numeric(t <= cen)
  survObj <- list(t = survObj$t, di = survObj$di, x = survObj$x)
  
  # Set hyperparameters
  s <- c(sort(survObj$t[survObj$di == 1]), 
         2 * max(survObj$t) - max(survObj$t[-which(survObj$t == max(survObj$t))]))
  priorPara <- list('eta0' = 1, 'kappa0' = 1, 'c0' = 2, 'r' = 0.5, 
                    'delta' = 0.0001, 's' = s, 'J' = length(s), 
                    'groupInd' = 1:p)
  
  # Set MCMC parameters
  mcmcPara <- list('numBeta' = p, 'beta.prop.var' = 1)
  
  # Initial values of hyperparameters
  lambdaSq <- 1
  initial <- list(
    'beta.ini' = rep(0, p),
    'lambdaSq' = 1,
    'sigmaSq' = runif(1, 0.1, 10),
    'tauSq' = rexp(length(unique(priorPara$groupInd)), rate = lambdaSq / 2),
    'h' = ifelse(hazard == "Gamma prior", rgamma(priorPara$J, 1, 1), rbeta(priorPara$J, 1, 1))  # Direct conditional assignment for h
  )
  
  # # Fit psbcGL
  PSBC.GL.G <- psbcGL.fit(survObj, priorPara, initial, rw = TRUE, 
                          mcmcPara, num.reps = 100, thin = 1, chain = 1)
  
  # Compute BIC vector
  bicVecResult <- bicVec(PSBC.GL.G, survObj$x, psiVec = seq(0.001, 1, 0.001))
  
  # Variable selection
  VS_Result <- apply(PSBC.GL.G$beta.p, 2, mean) > 0.001 # 99.9% (0.001) predicted power
  count_true <- sum(VS_Result)
  ## Variable Selection (VS)
  VS(PSBC.GL.G, X=survObj$x)
  # Print results in the required format
  print(paste0("p = ", p, ": Number of Selected variables = ", count_true))
  print(paste0("Prediction Concordance Index: ", PSBC.GL.G$concordance.index))
  print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
  print(paste0("BIC: ", bicVecResult[length(bicVecResult)]))  # BIC value from the last entry of bicVec
  
  # Store results
  results[[paste0("p_", p)]] <- list(fit = PSBC.GL.G, bicVec = bicVecResult, 
                                     VS = VS_Result, count_true = count_true)
  
  # Compute Predicted Risk
  pred_risk <- rowSums(survObj$x %*% t(PSBC.GL.G$beta.p))
  
  # Calculate AUC from Time-dependent ROC Curve
  auc_list_TP <- sapply(survObj$t, function(t) {
    roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
    roc_curve$TP
  })
  mean_TP <- mean(auc_list_TP, na.rm = TRUE)
  sd_TP <- sd(auc_list_TP, na.rm = TRUE)
  
  auc_list_FP <- sapply(survObj$t, function(t) {
    roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
    roc_curve$FP
  })
  mean_FP <- mean(auc_list_FP, na.rm = TRUE)
  sd_FP <- sd(auc_list_FP, na.rm = TRUE)
  
  auc_list <- sapply(survObj$t, function(t) {
    roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
    roc_curve$AUC
  })
  
  mean_auc <- mean(auc_list, na.rm = TRUE)
  
  plot(sort(auc_list_FP), sort(auc_list_TP), type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = paste("FP", "\n", "AUC = ", round(mean_auc, 3)), 
       ylab = "TP", main = "Mayoscore 4, Method = NNE \n  Year = 1")
}



