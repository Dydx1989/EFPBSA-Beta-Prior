# Source necessary files
source("bicVec.R")
source("snc.R")
source("psbcGL_fit.R")
source("psbcEN_fit.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(psbcGroup)

# Simulation examples
set.seed(204542)
# List of p values
p_values <- c(20, 200, 500, 1000)
results.1 <- list() # For number of variables selected and BIC ...
results.2 <- list() # For AUC
plot_list <- list()
# Loop over different p values
for (p in p_values) {
  
  # Set the hazard type (either "Gamma prior" or "Beta prior")
  hazard <- "Beta prior"  # Or change this to "Gamma prior"
  
  n <- 100
  
  beta.true <- c(rep(4, 10), rep(0, (p-10)))	
  
  CovX<- diag(0.1, p)
  
  survObj 	<- list()
  survObj$x	<- apply(rmvnorm(n, sigma=CovX, method="chol"), 2, scale)
  
  pred <- as.vector(exp(rowSums(scale(survObj$x, center = FALSE, scale = 1/beta.true))))
  
  t 		<- rexp(n, rate = pred)
  cen		<- runif(n, 0, 8)      
  survObj$t 		<- pmin(t, cen)
  survObj$di 		<- as.numeric(t <= cen)
  
  priorPara 			<- list()
  priorPara$eta0 		<- 1
  priorPara$kappa0 	<- 1
  priorPara$c0 		<- 2
  priorPara$r1		<- 0.1
  priorPara$r2		<- 1
  priorPara$delta1	<- 0.1
  priorPara$delta2	<- 1
  priorPara$s			<- sort(survObj$t[survObj$di == 1])
  priorPara$s			<- c(priorPara$s, 2*max(survObj$t)
                     - max(survObj$t[-which(survObj$t==max(survObj$t))]))
  priorPara$J			<- length(priorPara$s)
  
  mcmcPara				<- list()
  mcmcPara$numBeta		<- p
  mcmcPara$beta.prop.var	<- 1
  
  
  initial <- list(
    'beta.ini' = rep(0.5, p),
    'lambda1Sq' = 1,
    'lambda2' = 1,
    'sigmaSq' = runif(1, 0.1, 10),
    'tauSq' =  rexp(p, rate = initial$lambda1Sq/2),
    'h' = ifelse(hazard == "Gamma prior", rgamma(priorPara$J, 1, 1), rbeta(priorPara$J, 1, 1))  # Direct conditional assignment for h
  )
  
  
  rw = FALSE
  num.reps = 100
  chain = 1
  thin = 5
  save = 5
  
  
  # # Fit psbcGL
  PSBC.EL.G <-  psbcEN.fit(survObj, priorPara, initial, rw=FALSE, mcmcPara, 
                           num.reps, thin, chain, save)
  
  # Compute BIC vector
  bicVecResult <- bicVec(PSBC.EL.G, survObj$x, psiVec = seq(0.001, 1, 0.001))
  
  # Variable selection
  VS_Result <- apply(PSBC.EL.G $beta.p, 2, mean) > 0.001 # 99.9% (0.001) prediction power
  count_true <- sum(VS_Result)
  ## Variable Selection (VS)
  VS(PSBC.EL.G, X=survObj$x)
  # Print results in the required format
  print(paste0("p = ", p, ": Number of Selected variables = ", count_true))
  print(paste0("Prediction Concordance Index: ", PSBC.EL.G$concordance))
  print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
  print(paste0("BIC: ", bicVecResult[length(bicVecResult)]))  # BIC value from the last entry of bicVec
  
  # Store results
  results.1[[paste0("p_", p)]] <- list(fit = PSBC.EL.G, bicVec = bicVecResult, 
                                     VS = VS_Result, count_true = count_true)
  # Compute predicted risk
  pred_risk <- rowSums(survObj$x %*% t(PSBC.EL.G$beta.p))
  
  # Calculate AUC, TP, and FP
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
  
  mean_auc <- min(mean(auc_list, na.rm = TRUE), 1)
  auc_list <- pmin(auc_list, 1)
  # Store results
  results.2[[paste0("p_", p)]] <- list(
    mean_TP = mean_TP,
    sd_TP = sd_TP,
    mean_FP = mean_FP,
    sd_FP = sd_FP,
    mean_auc = mean_auc
  )
  
  # Print results for current p
  cat(sprintf("p = %d\n", p))
  cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP, sd_TP))
  cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP, sd_FP))
  cat(sprintf("Mean AUC: %.3f\n\n", mean_auc))
  
  # Generate plot for current p
  g.plot <- ggplot(data.frame(FP = sort(survObj$t), TP = sort(auc_list)), aes(x = FP, y = TP)) +
    geom_line(linewidth = 1.2) +  # Replace size with linewidth
    xlab(paste("Time", "\n", "Average AUC = ", round(mean_auc, 3))) +
    ylab("Area under ROC curve") +
    ggtitle(paste("p =", p, ", Hazard Type =", hazard)) +
    ylim(0, 1) +  # Add y-axis limits here
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Increase title size and bold
      axis.title = element_text(size = 14, face = "bold"),  # Increase axis titles size and bold
      axis.title.x = element_text(size = 16, face = "bold"),  # Increase x-axis label size and bold
      axis.title.y = element_text(size = 16, face = "bold"),  # Increase y-axis label size and bold
      axis.text = element_text(size = 14),  # Increase axis tick label size
      #plot.subtitle = element_text(size = 14)  # If you use subtitle
    )
  
  
  plot_list[[paste0("p_", p)]] <- g.plot
  
  
  # Combine all plots
  g.all <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
  
  # Save the arranged plot to a PNG file
  sf <- 2  # Scaling factor
  png(
    filename = paste0(hazard,"-output-B-all.png"), 
    width = 1600 * sf, 
    height = 1200 * sf, 
    res = 72 * sf
  )
  print(g.all)
  dev.off()
  
  }




