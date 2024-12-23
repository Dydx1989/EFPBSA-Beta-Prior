
source("setting.interval.R")
source("bic4.R")
source("loglh2.R")
source("snc.R")
source("UpdateRPrw.R")
source("UpdateBH.R")
source("UpdateTau.GL.R")
source("UpdateSigma.GL.R")
source("UpdateLambda.GL.R")
source("UpdateRP.R")
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("survivalROC")
library(survcomp)
library(survivalROC)

psbcGL.fit <- function(survObj, priorPara, initial, rw = FALSE, mcmcPara, num.reps, thin, chain = 1, save = 1000) {
  
  survObj$n <- n <- length(survObj$t)
  survObj$p <- p <- dim(survObj$x)[2]
  
  # Setting up prior parameters
  eta0 <- priorPara$eta0
  kappa0 <- priorPara$kappa0
  c0 <- priorPara$c0
  r <- priorPara$r
  delta <- priorPara$delta
  s <- priorPara$s
  J <- priorPara$J <- length(priorPara$s)
  groupInd <- priorPara$groupInd
  groupNo <- priorPara$groupNo <- unique(priorPara$groupInd)
  K <- priorPara$K <- length(groupNo)
  m_k <- priorPara$m_k
  
  m_k <- rep(NA, K)
  for (i in 1:K) {
    m_k[i] <- sum(groupInd == groupNo[i])
  }
  priorPara$m_k <- m_k
  
  intv <- setting.interval(survObj$t, survObj$di, priorPara$s, priorPara$J)
  priorPara$ind.r <- intv$ind.r
  priorPara$ind.d <- intv$ind.d
  priorPara$ind.r_d <- intv$ind.r_d
  priorPara$d <- intv$d
  
  ini <- initial
  beta.ini <- ini$beta.ini
  lambdaSq <- ini$lambdaSq
  sigmaSq <- ini$sigmaSq
  tauSq <- ini$tauSq
  h <- ini$h
  
  mcmcPara$beta.prop.me <- beta.ini
  
  tauSq.exp <- rep(NA, p)
  for (i in 1:K) {
    tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
  }
  ini$sd.be <- sqrt(sigmaSq * tauSq.exp)
  ini$xbeta <- as.vector(survObj$x %*% beta.ini)
  
  be.normSq <- c()
  for (i in 1:K) {
    be.normSq[i] <- sum(beta.ini[which(groupInd == groupNo[i])]^2)
  }
  ini$be.normSq <- be.normSq
  
  H.star <- alpha0 <- c()
  for (j in 1:J) {
    H.star[j] <- eta0 * s[j]^kappa0
    alpha0[j] <- c0 * H.star[j]
  }
  priorPara$hPriorSh <- diff(c(0, alpha0))
  
  ## For posterior samples
  mcmcOutcome <- list()
  mcmcOutcome$initial <- initial
  mcmcOutcome$priorPara <- priorPara
  
  beta.p <- beta.ini
  h.p <- h
  tauSq.p <- tauSq
  mcmcOutcome$sigmaSq.p <- sigmaSq
  mcmcOutcome$lambdaSq.p <- lambdaSq
  mcmcOutcome$accept.beta <- c(rep(0, p))
  
  outcomeSum <- list()
  dir.create('mcmcOutcome', showWarnings = FALSE)
  
  # MCMC Sampling
  for (M in 1:num.reps) {
    if (M %% 1000 == 0) {
      cat("Chain", chain, "Iteration", M, fill = TRUE)
    }
    
    # Update regression parameters
    sampleRP <- if (rw == FALSE) UpdateRP(survObj, priorPara, mcmcPara, ini) else UpdateRPrw(survObj, priorPara, mcmcPara, ini)
    beta.ini <- ini$beta.ini <- sampleRP$beta.ini
    xbeta <- ini$xbeta <- sampleRP$xbeta
    mcmcOutcome$accept.beta <- mcmcOutcome$accept.beta + sampleRP$accept
    
    for (i in 1:K) {
      be.normSq <- ini$be.normSq[i] <- sum(beta.ini[which(groupInd == groupNo[i])]^2)
    }
    
    h <- ini$h <- UpdateBH(survObj, priorPara, ini)
    tauSq <- ini$tauSq <- UpdateTau.GL(survObj, priorPara, ini)
    
    for (i in 1:K) {
      tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
    }
    
    sigmaSq <- ini$sigmaSq <- UpdateSigma.GL(survObj, priorPara, ini)
    lambdaSq <- ini$lambdaSq <- UpdateLambda.GL(survObj, priorPara, ini)
    ini$sd.be <- sqrt(sigmaSq * tauSq.exp)
    
    # Storing posterior samples
    if (M %% thin == 0) {
      beta.p <- rbind(beta.p, beta.ini, deparse.level = 0)
      h.p <- rbind(h.p, h, deparse.level = 0)
      tauSq.p <- rbind(tauSq.p, tauSq, deparse.level = 0)
      mcmcOutcome$sigmaSq.p <- c(mcmcOutcome$sigmaSq.p, sigmaSq)
      mcmcOutcome$lambdaSq.p <- c(mcmcOutcome$lambdaSq.p, lambdaSq)
      mcmcOutcome$ini <- ini
    }
    
    # Save MCMC outcomes
    if (M %% save == 0 || M == num.reps) {
      save(mcmcOutcome, file = paste("mcmcOutcome/otherAll.ch", chain, ".Rdata", sep = ""))
      save(beta.p, file = paste("mcmcOutcome/betaAll.ch", chain, ".Rdata", sep = ""))
      save(tauSq.p, file = paste("mcmcOutcome/tauSqAll.ch", chain, ".Rdata", sep = ""))
      save(h.p, file = paste("mcmcOutcome/hAll.ch", chain, ".Rdata", sep = ""))
    }
  }
  
  # Compute Concordance Index for Prediction Accuracy
  pred_risk <- rowSums(survObj$x %*% t(beta.p)) # posterior mean predicted risk for each individual, collapsing over all MCMC samples
  true_risk <- survObj$t
  concordance <- concordance.index(pred_risk, survObj$t, survObj$di)$c.index
  
  ret <- list(beta.p = beta.p, h.p = h.p, tauSq.p = tauSq.p, mcmcOutcome = mcmcOutcome, t = survObj$t, di = survObj$di, concordance = concordance)
  class(ret) <- "psbcGL"
  
  cat("Prediction Concordance Index:", concordance, "\n")
  
  return(ret)
}



# Fit psbcGL
#FitGL<- psbcGL.fit(survObj, priorPara, initial, rw = TRUE, 
  #             mcmcPara, num.reps = 100, thin = 1, chain = 1)

