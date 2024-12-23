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
DutchBC=read.csv("DBCD_Data.csv")
dim(DutchBC)
#data("nki70")
x<-DutchBC[,-c(1,2)]
p<-ncol(x)
#set.seed(123)
n = nrow(x)

###############################
y<-DutchBC[,c(1,2)]
t=y$time
di<-y$status
#set.seed(204542)
survObj_3 = list(t = t, di = di, x = scale(as.matrix(DutchBC[,-c(1,2)])))
p = ncol(survObj$x)

s = c(sort(survObj_3$t[survObj_3$di == 1]), 2 * max(survObj_3$t) - max(survObj_3$t[-which(survObj_3$t==max(survObj_3$t))]))
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
PSBC.GL.G.DutchBC = psbcGL.fit(survObj_3, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
#PSBC.GL.G.DutchBC$beta.p
save(PSBC.GL.G.DutchBC, file = "PSBC.GL.G.DutchBC.RData")
# Compute BIC vector
bicVecResult.GL.G <- bicVec(PSBC.GL.G.DutchBC, survObj_3$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result.GL.G <- apply(PSBC.GL.G.DutchBC$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true.GL.G <- sum(VS_Result.GL.G)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true.GL.G))
print(paste0("Prediction Concordance Index: ", PSBC.GL.G.DutchBC$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.GL.G[length(bicVecResult.GL.G)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj_3$x %*% t(PSBC.GL.G.DutchBC$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP.GL.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP.GL.G <- mean(auc_list_TP.GL.G, na.rm = TRUE)
sd_TP.GL.G <- sd(auc_list_TP.GL.G, na.rm = TRUE)

auc_list_FP.GL.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP.GL.G <- mean(auc_list_FP.GL.G, na.rm = TRUE)
sd_FP.GL.G <- sd(auc_list_FP.GL.G, na.rm = TRUE)

auc_list.GL.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

#save(auc_list.GL.G, file = "auc_list.GL.G.dutch.RData")
mean_auc.GL.G <- min(mean(auc_list.GL.G, na.rm = TRUE), 1)
auc_list.GL.G <- pmin(auc_list.GL.G, 1)
saveRDS(auc_list.GL.G ,"auc_list.GL.G.dutch.rds")

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP.GL.G, sd_TP.GL.G))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP.GL.G, sd_FP.GL.G))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc.GL.G))


### Beta prior 
#PSBC-GL-B
initial = list('beta.ini'= rep(0, p), 'lambdaSq' = 1, 'sigmaSq' = runif(1, 0.1, 10),
               'tauSq' = rexp(length(unique(priorPara$groupInd)), 'rate' = lambdaSq / 2),
               'h' = rbeta(priorPara$J, 1, 1))
# in real applications, 'num.reps' should be large enough (e.g. 20000, 40000) and 'chain' to be > 1
PSBC.GL.B.DutchBC = psbcGL.fit(survObj_3, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
#PSBC.GL.G.nki70$beta.p

# Compute BIC vector
bicVecResult.GL.B <- bicVec(PSBC.GL.B.DutchBC, survObj_3$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result.GL.B <- apply(PSBC.GL.B.DutchBC$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true.GL.B <- sum(VS_Result.GL.B)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true.GL.B))
print(paste0("Prediction Concordance Index: ", PSBC.GL.B.DutchBC$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.GL.B[length(bicVecResult.GL.B)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj_3$x %*% t(PSBC.GL.B.DutchBC$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP.GL.B <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP.GL.B <- mean(auc_list_TP.GL.B, na.rm = TRUE)
sd_TP.GL.B <- sd(auc_list_TP.GL.B, na.rm = TRUE)

auc_list_FP.GL.B <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP.GL.B <- mean(auc_list_FP.GL.B, na.rm = TRUE)
sd_FP.GL.B <- sd(auc_list_FP.GL.B, na.rm = TRUE)

auc_list.GL.B <- sapply(survObj$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc.GL.B <- min(mean(auc_list.GL.B, na.rm = TRUE), 1)
auc_list.GL.B <- pmin(auc_list.GL.B, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP.GL.B, sd_TP.GL.B))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP.GL.B, sd_FP.GL.B))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc.GL.B))



############# ENASTIC NET ###########
## Gamma prior 
#priorPara 			<- list()
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
PSBC.EN.G.DutchBC = psbcEN.fit(survObj_3, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)


# Compute BIC vector
bicVecResult.EN.G <- bicVec(PSBC.EN.G.DutchBC, survObj_3$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result.EN.G <- apply(PSBC.EN.G.DutchBC$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true.EN.G <- sum(VS_Result.EN.G)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true.EN.G))
print(paste0("Prediction Concordance Index: ", PSBC.EN.G.nki70$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.EN.G[length(bicVecResult.EN.G)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj_3$x %*% t(PSBC.EN.G.DutchBC$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP.EN.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP.EN.G <- mean(auc_list_TP.EN.G, na.rm = TRUE)
sd_TP.EN.G <- sd(auc_list_TP.EN.G, na.rm = TRUE)

auc_list_FP.EN.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP.EN.G <- mean(auc_list_FP.EN.G, na.rm = TRUE)
sd_FP.EN.G<- sd(auc_list_FP.EN.G, na.rm = TRUE)

auc_list.EN.G <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc.EN.G <- min(mean(auc_list.EN.G, na.rm = TRUE), 1)
auc_list.EN.G <- pmin(auc_list.EN.G, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP.EN.G, sd_TP.EN.G))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP.EN.G, sd_FP.EN.G))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc.EN.G))


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
PSBC.EN.B.DutchBC = psbcEN.fit(survObj_3, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)


# Compute BIC vector
bicVecResult.EN.B<- bicVec(PSBC.EN.B.DutchBC, survObj_3$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result.EN.B <- apply(PSBC.EN.B.nki70$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true.EN.B <- sum(VS_Result.EN.B)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true.EN.B))
print(paste0("Prediction Concordance Index: ", PSBC.EN.B.DutchBC$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.EN.B[length(bicVecResult.EN.B)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj_3$x %*% t(PSBC.EN.B.DutchBC$beta.p))

# Calculate AUC, TP, and FP
auc_list_TP.EN.B <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$TP
})
mean_TP.EN.B <- mean(auc_list_TP.EN.B, na.rm = TRUE)
sd_TP.EN.B <- sd(auc_list_TP.EN.B, na.rm = TRUE)

auc_list_FP.EN.B <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$FP
})
mean_FP.EN.B <- mean(auc_list_FP.EN.B, na.rm = TRUE)
sd_FP.EN.B <- sd(auc_list_FP.EN.B, na.rm = TRUE)

auc_list.EN.B <- sapply(survObj_3$t, function(t) {
  roc_curve <- survivalROC(Stime = survObj_3$t, status = survObj_3$di, marker = pred_risk, predict.time = t, method = "KM")
  roc_curve$AUC
})

mean_auc.EN.B <- min(mean(auc_list.EN.B, na.rm = TRUE), 1)
auc_list.EN.B <- pmin(auc_list.EN.B, 1)

cat(sprintf("Mean TP: %.3f, SD TP: %.3f\n", mean_TP.EN.B, sd_TP.EN.B))
cat(sprintf("Mean FP: %.3f, SD FP: %.3f\n", mean_FP.EN.B, sd_FP.EN.B))
cat(sprintf("Mean AUC: %.3f\n\n", mean_auc.EN.B))





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



##################################

auc_list.GL.G=readRDS("auc_list.GL.G.dutch.rds")
auc_list.GL.B=readRDS("auc_list.GL.B.dutch.rds")
auc_list.EN.G= readRDS("auc_list.EN.G.dutch.rds")
auc_list.EN.B= readRDS("auc_list.EN.B.dutch.rds")

Suv.Time=survObj_3$t
sorted_auc_list.GL.G <- auc_list.GL.G[order(auc_list.GL.G, na.last = TRUE)]
sorted_auc_list.GL.B <- auc_list.GL.B[order(auc_list.GL.B, na.last = TRUE)]
sorted_auc_list.EN.G <- auc_list.EN.G[order(auc_list.EN.G, na.last = TRUE)]
sorted_auc_list.EN.B <- auc_list.EN.B[order(auc_list.EN.B, na.last = TRUE)]

# Load necessary library
library(ggplot2)

# Combine the data into a data frame
dataPSBC <- data.frame(
  x = rep(Suv.Time, 4),
  y = c(sorted_auc_list.GL.G ,sorted_auc_list.GL.B,sorted_auc_list.EN.G,sorted_auc_list.EN.B ),
  line = rep(c("PSBC-GL-G","PSBC-GL-B", "PSBC-EN-G","PSBC-EN-B"), each = length(sorted_auc_list.GL.B))
)

# Create the ggplot
Plotdutch <- ggplot(dataPSBC, aes(x = x, y = y, color = line)) +
  geom_line(size = 1.5) +  # Increase line thickness
  scale_color_manual(values = c("blue", "red", "green", "black")) +
  labs(
    x = "Time",
    y = "Area under ROC curve (AUC)",
    color = "Methods"
  ) +
  theme_minimal(base_size = 15) +  # Increase base font size for theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Larger, bold title
    axis.title.x = element_text(size = 18),  # Larger x-axis title font
    axis.title.y = element_text(size = 18),  # Larger y-axis title font
    axis.text = element_text(size = 15),     # Larger axis labels
    legend.title = element_text(size = 16),  # Larger legend title
    legend.text = element_text(size = 14)    # Larger legend text
  )

# Display the plot
print(Plotdutch)

# Save the arranged plot to a PNG file
sf <- 2  # Scaling factor
datatype<- "Dutch-AUC"
png(
  filename = paste0(datatype,"-Plot.png"), 
  width = 700* sf, 
  height = 500 * sf, 
  res = 72 * sf
)
print(Plotdutch)
dev.off()

### nki70'



##################################

auc_list.GL.G=readRDS("auc_list.GL.G.nki70.rds")
auc_list.GL.B=readRDS("auc_list.GL.B.nki70.rds")
auc_list.EN.G= readRDS("auc_list.EN.G.nki70.rds")
auc_list.EN.B= readRDS("auc_list.EN.B.nki70.rds")

Suv.Time=survObj$t
sorted_auc_list.GL.G <- auc_list.GL.G[order(auc_list.GL.G, na.last = TRUE)]
sorted_auc_list.GL.B <- auc_list.GL.B[order(auc_list.GL.B, na.last = TRUE)]
sorted_auc_list.EN.G <- auc_list.EN.G[order(auc_list.EN.G, na.last = TRUE)]
sorted_auc_list.EN.B <- auc_list.EN.B[order(auc_list.EN.B, na.last = TRUE)]

# Load necessary library
library(ggplot2)

# Combine the data into a data frame
dataPSBC <- data.frame(
  x = rep(Suv.Time, 4),
  y = c(sorted_auc_list.GL.G ,sorted_auc_list.GL.B,sorted_auc_list.EN.G,sorted_auc_list.EN.B ),
  line = rep(c("PSBC-GL-G","PSBC-GL-B", "PSBC-EN-G","PSBC-EN-B"), each = length(sorted_auc_list.GL.B))
)

# Create the ggplot
PlotNKI70 <- ggplot(dataPSBC, aes(x = x, y = y, color = line)) +
  geom_line(size = 1.5) +  # Increase line thickness
  scale_color_manual(values = c("blue", "red", "green", "black")) +
  labs(
    x = "Time",
    y = "Area under ROC curve (AUC)",
    color = "Methods"
  ) +
  theme_minimal(base_size = 15) +  # Increase base font size for theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Larger, bold title
    axis.title.x = element_text(size = 18),  # Larger x-axis title font
    axis.title.y = element_text(size = 18),  # Larger y-axis title font
    axis.text = element_text(size = 15),     # Larger axis labels
    legend.title = element_text(size = 16),  # Larger legend title
    legend.text = element_text(size = 14)    # Larger legend text
  )

# Display the plot
print(PlotNKI70)


all.plot=ggarrange(PlotNKI70,                                                 # First row with scatter plot
          ggarrange(Plotdutch, Plotchop, ncol = 2, labels = c("B: DutchBC data", "C: DLBCL data")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A: NKI70 data"                                        # Labels of the scatter plot
) 

# Save the arranged plot to a PNG file
sf <- 2  # Scaling factor
datatype<- "All.data-AUC"
png(
  filename = paste0(datatype,"-Plot.png"), 
  width = 1000* sf, 
  height = 700 * sf, 
  res = 72 * sf
)
print(all.plot)
dev.off()

