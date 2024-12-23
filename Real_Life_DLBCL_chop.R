library("dplyr")
library("ggplot2")
library("survival")
library("survminer")
#library("M3C")
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
# Source necessary files
source("bicVec.R")
source("snc.R")
source("psbcGL_fit.R")

# "PSBC-GL-G", "PSBC-EN-G", "PSBC-GL-B", "PSBC-EN-B

Datachop=read.csv("chop.csv")
Datachop=Datachop[,-1]
dim(Datachop)
head(Datachop[,1:5])
set.seed(123)
survObj = list(t =Datachop[,1] , di = Datachop[,2], x = as.matrix(Datachop[,-c(1,2)]))
p = ncol(survObj$x)
# set hyperparameters. 
# For Lasso prior (i.e. 'groupInd'= 1:p), larger ratio r/delta tends to force the posterior betas to be more concentrated at 0
# For group Lasso prior (i.e. 'groupInd' as group indicator for covariates), larger ratio r/delta tends to force stronger grouping effect of covariates
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

PSBC.GL.G = psbcGL.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
save(PSBC.GL.G, file = "PSBC.GL.G.DLBCL.RData")
# Compute BIC vector
bicVecResult <- bicVec(PSBC.GL.G, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result <- apply(PSBC.GL.G$beta.p, 2, mean) > mean(PSBC.GL.G$beta.p) # 99.9% (0.001) prediction power
count_true <- sum(VS_Result)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true))
print(paste0("Prediction Concordance Index: ", PSBC.GL.G$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult[length(bicVecResult)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.GL.G$beta.p))

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
PSBC.GL.B= psbcGL.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
save(PSBC.GL.B, file = "PSBC.GL.B.DLBCL.RData")

#PSBC.GL.G.nki70$beta.p

# Compute BIC vector
bicVecResult.GL.B <- bicVec(PSBC.GL.B, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_2 <- apply(PSBC.GL.B$beta.p, 2, mean) > 0.02 # 99.9% (0.001) prediction power
count_true_2 <- sum(VS_Result_2)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_2))
print(paste0("Prediction Concordance Index: ", PSBC.GL.B$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult.GL.B[length(bicVecResult.GL.B)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.GL.B$beta.p))

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
PSBC.EN.G = psbcEN.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
save(PSBC.EN.G, file = "PSBC.EN.G.DLBCL.RData")

# Compute BIC vector
bicVecResult_3 <- bicVec(PSBC.EN.G, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_3 <- apply(PSBC.EN.G$beta.p, 2, mean) > 0.02 # 99.8% (0.001) prediction power
count_true_3 <- sum(VS_Result_3)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_3))
print(paste0("Prediction Concordance Index: ", PSBC.EN.G$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult_3[length(bicVecResult_3)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.EN.G$beta.p))

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
PSBC.EN.B = psbcEN.fit(survObj, priorPara, initial, rw = TRUE, mcmcPara, num.reps = 100, thin = 1, chain = 1)
save(PSBC.EN.B, file = "PSBC.EN.B.DLBCL.RData")

# Compute BIC vector
bicVecResult_4 <- bicVec(PSBC.EN.B, survObj$x, psiVec = seq(0.001, 1, 0.001))

# Variable selection
VS_Result_4 <- apply(PSBC.EN.B$beta.p, 2, mean) > 0.2 # 99.9% (0.001) prediction power
count_true_4 <- sum(VS_Result_4)
## Variable Selection (VS)
#VS(PSBC.GL.G, X=survObj$x)
# Print results in the required format
print(paste0("p = ", p, ": Number of Selected variables = ", count_true_4))
print(paste0("Prediction Concordance Index: ", PSBC.EN.B$concordance))
#print("SNC: 1000 out of 1000")  # Assuming SNC is always 1000 out of 1000 based on the output
print(paste0("BIC: ", bicVecResult_4[length(bicVecResult_4)]))  # BIC value from the last entry of bicVec

# Compute predicted risk
pred_risk <- rowSums(survObj$x %*% t(PSBC.EN.B$beta.p))

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












# burn-in the first half MCMC iterations
beta_p = BayesLassofit$beta.p[-(1:51), ]
beta_mean = colMeans(beta_p)
beta_L = apply(beta_p, 2, quantile, 0.025)
beta_U = apply(beta_p, 2, quantile, 0.975)
tbl = data.frame(term =colnames(Datachop[,-c(1,2)]), estimate = beta_mean,  conf.low = beta_L,  conf.high = beta_U)
tbl$term = factor(tbl$term, levels = tbl$term)
## ESTIMATE, CREDIBLE REGION
tbl
GGally::ggcoef(tbl) + xlab(expression(Posterior~~beta)) + ylab("")


Datachop=read.csv("chop.csv")
Datachop=Datachop[,-1]
dim(Datachop)
head(Datachop[,1:5])
set.seed(123)
survObj = list(t =Datachop[,1] , di = Datachop[,2], x = as.matrix(Datachop[,-c(1,2)]))
## for demonstration simplicity, PAM50 genes are used here
x = as.matrix(Datachop[,-c(1,2)])
Datachop[,1][Datachop[,1] == 0] <- 1
y = cbind(time=Datachop[,1], status=Datachop[,2])
set.seed(123)
n = nrow(survObj$x)
idx = sample(1:n, n * 0.8, replace = FALSE)
x_train = survObj$x[idx, ]
y_train = y[idx, ]

x_validate = x[-idx, ]
y_validate = y[-idx, ]

# train a Lasso Cox model, similarly for other Cox-type models
set.seed(123)
library(caret)
# set penalty factor without penalizing the two demographical variables
pf = c(rep(0, 2), rep(1, ncol(x) - 2))
cvfit = glmnet::cv.glmnet(x_train, y_train, family = "cox", nfolds = 5, penalty.factor = pf)
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



####################################
auc_list
auc_list_2
auc_list_3
auc_list_4 

auc_list.GL.G=auc_list
auc_list.GL.B=auc_list_2
auc_list.EN.G= auc_list_3
auc_list.EN.B= auc_list_4 

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
Plotchop <- ggplot(dataPSBC, aes(x = x, y = y, color = line)) +
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
print(Plotchop)

# Save the arranged plot to a PNG file
sf <- 2  # Scaling factor
datatype<- "DLBCL-AUC"
png(
  filename = paste0(datatype,"-Plot.png"), 
  width = 800* sf, 
  height = 600 * sf, 
  res = 72 * sf
)
print(Plotchop)
dev.off()
