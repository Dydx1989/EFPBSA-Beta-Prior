#library(remotes)
#remotes::install_github("nyiuab/BhGLM")
library("dplyr")
library("ggplot2")
library("survival")
library("survminer")
library("M3C")
library("glmnet")
library("plotmo")
library("grpreg")
library("SGL")
library("psbcGroup")
library("GGally")
library("BhGLM")
library("risksetROC")
library("riskRegression")
library("peperr")
library("c060")
library("rms")
library("survAUC")
library("regplot")
############
# Source necessary files
source("bicVec.R")
source("snc.R")
source("psbcGL_fit.R")
source("psbcEN_fit.R")
#####################################################
#library(penalized)
nki70=read.csv("nki70_genes.csv")
dim(nki70)
head(nki70)
#data("nki70")
x<-nki70[,-c(1,2)]
p<-ncol(x)
#set.seed(123)
n = nrow(x)
dim(nki70)
###############################
y<-nki70[,c(1,2)]
t=y$time
di<-y$event
set.seed(204542)
p=75
n=144
survObj = list(t = t, di = di, x = as.matrix(nki70[,-c(1,2)]))
p = ncol(survObj$x)

s = c(sort(survObj$t[survObj$di == 1]), 2 * max(survObj$t) - max(survObj$t[-which(survObj$t==max(survObj$t))]))
priorPara = list('eta0' = 1, 'kappa0' = 1, 'c0'= 2, 'r' = 0.5, 
                 'delta' = 0.0001, 's'= s, 'J'=length(s), 'groupInd'= 1:p)
# set MCMC parameters
mcmcPara = list('numBeta'= p, 'beta.prop.var'= 1)
# set initial values of hyperparameters
lambdaSq = 1
initial = list('beta.ini'= rep(0, p), 'lambdaSq' = 1, 'sigmaSq' = runif(1, 0.1, 10),
               'tauSq' = rexp(length(unique(priorPara$groupInd)), 'rate' = lambdaSq / 2),
               'h' = rgamma(priorPara$J, 1, 1))
# in real applications, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
PSBC.GL.G.nki70 = psbcGL.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
PSBC.GL.G.nki70$beta.p

# Compute BIC vector
bicVecResult <- bicVec(PSBC.GL.G.nki70, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result <- apply(PSBC.GL.G.nki70$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true <- sum(VS_Result)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true))
print(paste0("Prediction Concordance Index: ", PSBC.GL.G.nki70$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult[length(bicVecResult)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.GL.G.nki70$beta.p))

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

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP, sd_TP))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP, sd_FP))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc))


### Beta prior 
#PSBC-GL-B
initial = list('beta.ini'= rep(0, p), 'lambdaSq' = 1, 'sigmaSq' = runif(1, 0.1, 10),
               'tauSq' = rexp(length(unique(priorPara$groupInd)), 'rate' = lambdaSq / 2),
               'h' = rbeta(priorPara$J, 1, 1))
# in real applications, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
PSBC.GL.B.nki70 = psbcGL.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
#PSBC.GL.G.nki70$beta.p

# Compute BIC vector
bicVecResult.GL.B <- bicVec(PSBC.GL.B.nki70, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_2 <- apply(PSBC.GL.B.nki70$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true_2 <- sum(VS_Result_2)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_2))
print(paste0("Prediction Concordance Index: ", PSBC.GL.B.nki70$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.GL.B[length(bicVecResult.GL.B)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.GL.B.nki70$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP_2 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP_2 <- mean(auc_list_TP_2, na.rm = TRUE)
sd_TP_2 <- sd(auc_list_TP_2, na.rm = TRUE)

auc_list_FP_2 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP_2 <- mean(auc_list_FP_2, na.rm = TRUE)
sd_FP_2 <- sd(auc_list_FP_2, na.rm = TRUE)

auc_list_2 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc_2 <- min(mean(auc_list_2, na.rm = TRUE), 1)
auc_list_2 <- pmin(auc_list_2, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP_2, sd_TP_2))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP_2, sd_FP_2))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc_2))



############# ENASTIC NET ###########
## Gamma prior 
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
  'h' = rgamma(priorPara$J, 1, 1)  # Direct conditional assignment for h
)

# in real applications, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
PSBC.EN.G.nki70 = psbcEN.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)


# Compute BIC vector
bicVecResult_3 <- bicVec(PSBC.EN.G.nki70, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_3 <- apply(PSBC.EN.G.nki70$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true_3 <- sum(VS_Result_3)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_3))
print(paste0("Prediction Concordance Index: ", PSBC.EN.G.nki70$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult_3[length(bicVecResult_3)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.EN.G.nki70$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP_3 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP_3 <- mean(auc_list_TP_3, na.rm = TRUE)
sd_TP_3 <- sd(auc_list_TP_3, na.rm = TRUE)

auc_list_FP_3 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP_3 <- mean(auc_list_FP_3, na.rm = TRUE)
sd_FP_3 <- sd(auc_list_FP_3, na.rm = TRUE)

auc_list_3 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc_3 <- min(mean(auc_list_3, na.rm = TRUE), 1)
auc_list_3 <- pmin(auc_list_3, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP_3, sd_TP_3))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP_3, sd_FP_3))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc_3))


### Beta Prior: Elastic net

initial <- list(
  'beta.ini' = rep(0.5, p),
  'lambda1Sq' = 1,
  'lambda2' = 1,
  'sigmaSq' = runif(1, 0.1, 10),
  'tauSq' =  rexp(p, rate = initial$lambda1Sq/2),
  'h' = rbeta(priorPara$J, 1, 1)  # Direct conditional assignment for h
)

# in real applications, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
PSBC.EN.B.nki70 = psbcEN.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)


# Compute BIC vector
bicVecResult_4 <- bicVec(PSBC.EN.B.nki70, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_4 <- apply(PSBC.EN.B.nki70$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true_4 <- sum(VS_Result_4)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_4))
print(paste0("Prediction Concordance Index: ", PSBC.EN.B.nki70$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult_4[length(bicVecResult_4)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.EN.B.nki70$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP_4 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP_4 <- mean(auc_list_TP_4, na.rm = TRUE)
sd_TP_4 <- sd(auc_list_TP_4, na.rm = TRUE)

auc_list_FP_4 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP_4 <- mean(auc_list_FP_4, na.rm = TRUE)
sd_FP_4 <- sd(auc_list_FP_4, na.rm = TRUE)

auc_list_4 <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj$t, status = survObj$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc_4 <- min(mean(auc_list_4, na.rm = TRUE), 1)
auc_list_4 <- pmin(auc_list_4, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP_4, sd_TP_4))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP_4, sd_FP_4))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc_4))





#t, x, alpha, beta, beta_all, sigmaSq, tauSq, nu0, sigSq0, alpha0, h0
#bicVec(PSBC.GL.G.nki70, X=survObj$x)
#Lambda<-PSBC.GL.G.nki70$mcmcOutcome$lambdaSq.p
#hist(Lambda, main = "p = 77, r = 0.5, #nonzero=2")

############################################
# burn-in the first half MCMC iterations
beta_p = BayesLassofit$beta.p[-(1:51), ]
beta_mean = colMeans(beta_p)
beta_L = apply(beta_p, 2, quantile, 0.025)
beta_U = apply(beta_p, 2, quantile, 0.975)
tbl = data.frame(term =colnames(survObj$x), estimate = beta_mean,  conf.low = beta_L,  conf.high = beta_U)
tbl$term = factor(tbl$term, levels = tbl$term)
## ESTIMATE, CREDIBLE REGION
tbl
GGally::ggcoef(tbl) + xlab(expression(Posterior~~beta)) + ylab("")


# Bayesian Cox model with Elastic Net prior
set.seed(123)
# set hyperparameters
# Larger ratio r1/delta1 forces the posterior betas to be more concentrated at 0
# Larger ratio r2/delta2 forces stronger grouping effect of covariates
priorPara = list('eta0' = 1, 'kappa0' = 1, 'c0'= 2, 'r1' = 0.1, 'r2' = 1, 
                 'delta1' = 0.1, 'delta2' = 1, 's'= s, 'J' = length(s))
# set MCMC parameters
mcmcPara = list('numBeta'= p, 'beta.prop.var'= 1)
# set initial values of hyperparameters
initial = list('beta.ini'= rep(0, p), 'lambda1Sq' = 1, 'lambda2' = 1, 'sigmaSq' = runif(1, 0.1, 10),
               'tauSq' = rexp(p, rate = 1 / 2), 'h' = rgamma(priorPara$J, 1, 1))
# in real application, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
BayesENfit = psbcEN(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
VS(BayesENfit, X=survObj$x)
#t, x, alpha, beta, beta_all, sigmaSq, tauSq, nu0, sigSq0, alpha0, h0
bicVec(BayesENfit, X=survObj$x)
Lambda<-BayesENfit$mcmcOutcome$lambda2.p
hist(Lambda, main = "p = 4922, r = 0.5, #nonzero=1")



# burn-in the first half MCMC iterations
EN_beta_p = BayesENfit$beta.p[52:101, ]
EN_beta_mean = colMeans(EN_beta_p)
EN_beta_L = apply(EN_beta_p, 2, quantile, 0.025)
EN_beta_U = apply(EN_beta_p, 2, quantile, 0.975)
EN_tbl = data.frame(term = colnames(x), estimate = EN_beta_mean, conf.low = EN_beta_L, conf.high = EN_beta_U)
EN_tbl$term = factor(EN_tbl$term, levels = EN_tbl$term)
## ESTIMATE, CREDIBLE REGION OR INTERVAL
EN_tbl
GGally::ggcoef(EN_tbl) + xlab(expression(Posterior~~beta)) + ylab("")


########### Model evaluation (classic) ###############
######################################################
######################################################
y<-Surv(nki70$time,nki70$event)
y<-Surv(DBCD$time,DBCD$status)
set.seed(123)
n = nrow(x)
xx<-as.matrix(x)
idx = sample(1:n, n * 0.8, replace = FALSE)
x_train = xx[idx, ]
y_train = y[idx, ]
x_validate = xx[-idx, ]
y_validate = y[-idx, ]
# train a Lasso Cox model, similarly for other Cox-type models
set.seed(123)
#class(x_train)
cvfit = cv.glmnet(x_train, y_train, family = "cox", nfolds = 5, penalty.factor = rep(2, ncol(x)))
pred_lp = predict(cvfit, newx = x_validate, s = cvfit$lambda.min, type = "link")

# dichotomize by prognostic scores (linear predictor)  by median to divide the validation patients into two groups
group_dichotomize = as.numeric(pred_lp > median(pred_lp))

# draw two survival curves based on KM estimation and compare them by a log-rank test
dat_tmp = data.frame(time = y_validate[, 1], status = y_validate[, 2], group = group_dichotomize)
sfit = survfit(Surv(time, status) ~ group, data = dat_tmp)

ggsurv = ggsurvplot(sfit, conf.int = TRUE, risk.table = TRUE, 
                    xlab = "Time since diagnosis (year)", legend = c(.2,.3),
                    legend.labs = c("Low risk", "High risk"), legend.title = "Dichotomized groups",  
                    risk.table.y.text.col = TRUE, risk.table.y.text = FALSE)
ggsurv$plot = ggsurv$plot + 
  annotate("text", x = 2.6, y = .03, label= paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))
ggsurv$table = ggsurv$table + labs(y = "Dichotomized\n groups")
ggsurv

group = pred_lp
group[pred_lp >= quantile(pred_lp, 2/3)] = 3
group[pred_lp >= quantile(pred_lp, 1/3) & pred_lp < quantile(pred_lp, 2/3)] = 2
group[pred_lp < quantile(pred_lp, 1/3)] = 1

# draw two survival curves based on KM estimation and compare them by a log-rank test
dat_tmp = data.frame(time = y_validate[, 1], status = y_validate[, 2], group = group)
sfit = survfit(Surv(time, status) ~ group, data = dat_tmp)

ggsurv = ggsurvplot(sfit, conf.int = TRUE, risk.table = TRUE, 
                    xlab = "Time since diagnosis (year)", legend = c(.2,.3),
                    legend.labs = c("Low risk", "Middle risk", "High risk"), legend.title = "Groups",  
                    risk.table.y.text.col = TRUE, risk.table.y.text = FALSE)
ggsurv$plot = ggsurv$plot + 
  annotate("text", x = 3.5, y = .05, label= paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))
ggsurv

###### ROC curve


ROC = risksetROC(Stime = y_validate[, 1], status = y_validate[, 2],
                 marker = pred_lp, predict.time = 5, method = "Cox", 
                 main = "ROC Curve", col = "seagreen3", type = "s", 
                 lwd = 2, xlab="1 - Specificity", ylab="Sensitivity") 
text(0.7, 0.2, paste("AUC =", round(ROC$AUC, 3)))



# Generate plot for current p
survObj_t= survObj$t
Auc.data_1=data.frame(FP = survObj$t, TP = auc_list)
Auc.data_2=data.frame(FP = survObj$t, TP = auc_list_2)
Auc.data_3=data.frame(FP = survObj$t, TP = auc_list_3)
Auc.data_4=data.frame(FP = survObj$t, TP = auc_list_4)
g.plot <- ggplot(Auc.data_1, aes(x = FP, y = TP)) +
  geom_line(linewidth = 1.2) +  # Replace size with linewidth
  xlab(paste("Time", "\n", "Average AUC = ", round(mean_auc, 3))) +
  ylab("Area under ROC curve") +
  #ggtitle(paste("p =", p, ", Hazard Type =", hazard)) +
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

