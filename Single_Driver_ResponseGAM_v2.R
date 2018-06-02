############
rm(list = ls())

set.seed(627)


############
# PACKAGES #
############
library(ggplot2)
library(meboot)
library(mgcv)
library(reshape2)
library(mvtnorm)
library(MASS)
library(nlme)
library(lme4)


######################
# Load functions to compare models 
######################

AICc <- function(mod) {
  K.c <- mod$rank
  N.c <- length(mod$residuals)
  AIC.c <- round(mod$aic + (2*K.c*(K.c+1)/(N.c-K.c-1)),3)
  return(AIC.c)
}



######################
# Make up data frame 
######################

tp_func <- function(x, a = 0.01, h = 10) {
  num <- a * x/(1 + a * x * h)
}

mydriver <- 0:100
myresponse <- sapply(mydriver, tp_func)
year <- 1950:(1950 + length(mydriver) - 1)
year <- sample(year, replace = FALSE)
ts.length <- year


myresponse <- myresponse + rnorm(length(myresponse), 0, 0.025)
plot(mydriver, myresponse)


df <- data.frame("year" = year ,"mydriver" = mydriver, "myresponse" = myresponse)

ggplot(df, aes(x=mydriver,y=myresponse))+
  geom_point()+
  #geom_smooth()+
  theme_classic()

ggplot(df, aes(x=year,y=myresponse))+
  geom_line()+
  #geom_smooth()+
  theme_classic()


#####################
# Fit the GAM model #
#####################

sp.len <- 200 # Spline length
nb <- 1000  # Number of bootstrap replicates
gam.mat <- matrix(nrow = sp.len, ncol = nb) ## create a matrix to be filled with bootstrapped splines (200 x 1000)
dif1 <- 1 # First derivative
dif2 <- 2 # Second derivative
ks <- 4   #If ks=3, no relationships have edf values > 2.0
rand <- rep(1, length(myresponse))

#Step 1. Fit GAM to answer whether temporal autocorrelation is important? Use the residuals from the gam and a log likelihood ratio test to calculate the "P.ac" value. A significant p.ac value suggests a model with autocorrelated error structure explains more of the variation in the residuals of the gam model than a model without autocorrelated error structure. Thus, using a GAMM is necessary to account for autocorrelation in the time series...use GAMM if p.ac < 0.05. If p.ac > 0.05, then GAM will be used to look for non-linearities, so this code will also fit the model and provide selection criteria (i.e. edf, GCV and AIC scores from GAM and Linear model (linear) to compare)

gam1  <- gam(myresponse ~ s(mydriver, bs = "tp", k = ks), optimMmethod = "GCV.Cp", se = T)
linear <- gam(myresponse ~ mydriver, method = "GCV.Cp", se = T)
dev.resid <- residuals(gam1, type = 'deviance')
lme1 <- lme(dev.resid ~ 1, random = ~1|rand, correlation = corAR1(form = ~year), method = 'ML')
lm1 <- lm(dev.resid ~ 1)
p.ac <- 1 - pchisq(2 * (logLik(lme1)[1] - logLik(lm1)[1]), 2)

#A negative value means the GAM with a smoother is a better model than the linear model  
delta.GCV.gam.lm <- summary(gam1)$sp.criterion - summary(linear)$sp.criterion
#A negative value means the GAM with a smoother is a better model than the linear model
delta.AIC.gam.lm <- AICc(gam1) - AICc(linear)
dev.diff.gam.lm <- summary(gam1)$dev.expl-summary(linear)$dev.expl

#Step 2. Fit GAMM to get selection criteria for relationships where p.ac < 0.05 (i.e. edf and AIC) and calculate deviance explained by GAMM ("gamm.dev.expl" below)
try(gamm1 <- gamm(myresponse ~ s(mydriver, bs= "tp",k = ks), optimMmethod="GCV.Cp",
                  se = T,correlation=corAR1(form=~ts.length)))
if (length(gamm1)==0) {
  print(imod)
  gamm1 <- gamm(myresponse ~ s(mydriver), se = T, correlation = corAR1(form = ~ts.length))
}

#Fit null model to compute deviance explained by gamm
null_model  <- glmmPQL(myresponse ~ 1, random = list(rand = ~1), family = 'gaussian')  
dr <- sum(residuals(gamm1$gam)^2)
dn0 <- sum(residuals(null_model)^2)
gamm.dev.expl <- (dn0-dr)/dn0


#Step 3. Fit linear model with autocorrelation (LMAC) to get selection criteria (i.e. AIC) and calculate deviance explained by LMAC ("lmac.dev.expl" below).
try(lmac <- gamm(myresponse ~ mydriver, optimMmethod = "GCV.Cp",
                 se = T, random =list(rand = ~1),
                 correlation =corAR1(form = ~ts.length)))
if (length(lmac)==0) {
  print(imod)
  lmac <- gamm(myresponse ~ mydriver, random =list(rand = ~1), 
               se = T,
               correlation = corAR1(form = ~ts.length))
}  

dr2 <- sum(residuals(lmac$gam)^2)
lmac.dev.expl <- (dn0-dr2)/dn0

# Calculate difference in deviance and AIC between GAMM and LMAC ("dev.diff.gamm.lmac" and "delta.AIC.gamm.lmac" respectively below).
dev.diff.gamm.lmac <- gamm.dev.expl-lmac.dev.expl
delta.AIC.gamm.lmac <- summary(gamm1$lme)$AIC-summary(lmac$lme)$AIC   #A negative value means the GAMM is a better model than the linear model with temporal autocorrelation



## Pull out relevant model outputs:
#FOR GAMM:
summary.gamm1 <- as.data.frame(cbind("GAMM",                    # Model name
                                     "myresponse",              # Response variable
                                     "mydriver",                # Pressure variable
                                     summary(gamm1$lme)$AIC,    # AICc
                                     summary(gamm1$lme)$logLik, # log likelihood
                                     gamm.dev.expl,             # Deviance explained by gamm
                                     summary(gamm1$gam)$edf,    # estimated degrees of freedom
                                     summary(gamm1$gam)$s.pv,   # p-value of the smoother
                                     summary(gamm1$gam)$r.sq,   # R-squared value
                                     NA,                        # placeholder for GCV value used in gam
                                     NA,                        # placeholder for delta.GCV data used as selection criteria for gams
                                     as.numeric(p.ac),          # p-value of log likelihood ratio test that tests whether the residuals from the gam are better explained by a model with or without temporal autocorrelation; p.ac values < 0.05 suggest autocorrelation is important and a gamm (instead of a gam) should be used.
                                     delta.AIC.gamm.lmac,       # delta AIC between gamm and lmac models; negative value means gamm is a better model than LMAC
                                     dev.diff.gamm.lmac         # difference in deviance explained between gamm and lmac
))

colnames(summary.gamm1) <- c("MODEL", "RESPONSE", "DRIVER", "AICc", "logLik","dev.expl", "edf", "Pvalue","R-squared","GCV","delta.GCV","p.ac","delta.AIC","diff.dev.expl")



#FOR Linear Model with Auto Correlation
summary.lmac <- as.data.frame(cbind("LMAC",                         # Model name
                                    "myresponse",                   # Response variable 
                                    "mydriver",                     # Pressure variable
                                    summary(lmac$lme)$AIC,          # AICc
                                    summary(lmac$lme)$logLik,       # log likelihood
                                    lmac.dev.expl,                  # Deviance explained by LMAC
                                    summary(lmac$gam)$residual.df,  # residual degrees of freedom
                                    summary(lmac$gam)$p.pv[[2]],    # p-value of the null hypothesis
                                    summary(lmac$gam)$r.sq,         # R-squared value
                                    NA,                             # placeholder for GCV value used in gam
                                    NA,                             # placeholder for delta.GCV data used as selection criteria for gams
                                    NA,                             # placeholder for p.ac for gamm
                                    NA,                             # placeholder for delta AIC                                      
                                    NA                              # placeholder for diff in deviance explained
))                          

colnames(summary.lmac)<- c("MODEL", "RESPONSE", "DRIVER", "AICc", "logLik","dev.expl", "edf", "Pvalue","R-squared","GCV","delta.GCV","p.ac","delta.AIC","diff.dev.expl")

#FOR GAM
summary.gam1 <- as.data.frame(cbind("GAM",                          # Model name
                                    "myresponse",                   # Response variable 
                                    "mydriver",                     # Pressure variable
                                    AICc(gam1),                     # AICc
                                    logLik(gam1),                   # Log likelihood
                                    summary(gam1)$dev.expl,         # deviance explained by gam
                                    summary(gam1)$edf,              # estimated degrees of freedom
                                    summary(gam1)$s.pv,             # p-value of the smoother
                                    summary(gam1)$r.sq,             # R-squared value
                                    summary(gam1)$sp.criterion,     # GCV value for gam
                                    delta.GCV.gam.lm,               # GAM GCV score minus LM GCV score; #A negative value means the GAM with a smoother is a better model than the linear model  
                                    NA,                             # placeholder for p.ac value for gamm
                                    delta.AIC.gam.lm,               # delta AICc between gam and linear model; a negaive value means the GAM with a smoother is a better model than the linear model
                                    dev.diff.gam.lm                 # Difference in deviance explained between the gam and the linear model
))
rownames(summary.gam1) <- NULL
colnames(summary.gam1)<- c("MODEL", "RESPONSE", "DRIVER", "AICc", "logLik","dev.expl", "edf", "Pvalue","R-squared","GCV","delta.GCV","p.ac","delta.AIC","diff.dev.expl")

#FOR LINEAR MODEL
summary.linear <- as.data.frame(cbind("Linear",                     # Model name
                                      "myresponse",                 # Response variable 
                                      "mydriver",                   # Pressure variable
                                      AICc(linear),                 # AICc
                                      logLik(linear),               # Log likelihood
                                      summary(linear)$dev.expl,     # deviance explained by linear model
                                      summary(linear)$residual.df,  # residual degrees of freedom
                                      summary(linear)$p.pv[[2]],    # p-value of the null hypothesis
                                      summary(linear)$r.sq,         # R-squared value
                                      summary(linear)$sp.criterion, # GCV value for linear model
                                      NA,                           # placeholder for delta.GCV data used as selection criteria for gams
                                      NA,                           # placeholder for p.ac for gamm
                                      NA,                           # placeholder for delta AIC
                                      NA                            # placeholder for diff in deviance explained
))                          

rownames(summary.linear) <- NULL
colnames(summary.linear) <- c("MODEL", "RESPONSE", "DRIVER", "AICc", "logLik","dev.expl", "edf", "Pvalue","R-squared","GCV","delta.GCV","p.ac","delta.AIC","diff.dev.expl")



### Identify "best model"#####
### 1a) Is p.ac <= 0.05? If yes, keep and evaluate selection criteria between GAMM and LMAC as best model.If no, move to step 2.
### 1b) Is GAMM edf > 2.0 for GAMM? If yes, keep. If no, LMAC is best.
### 1c) Is delta.AIC > 2.0 between GAMM and LMAC. If yes, GAMM is best. If no, LMAC is most parsimonious.
### 2a) If p.ac > 0.05, then revert to GAM model and ask if edf of GAM > 2.0? If yes, keep GAM. If no, linear model is best model.
### 2b) Is GCV minimized in GAM compared to Linear model? If delta.GCV.gam.lm is negative then keep GAM. If delta.GCV.gam.lm is positive then linear model is best.
### 2c) Is deltaAIC > 2.0 for GAM? If yes, then GAM is best model. If no, then linear model is best model.

summary.gamm1$best.model = ifelse(as.numeric(as.character(summary.gamm1$p.ac)) <= 0.05,
                                  ifelse(as.numeric(as.character(summary.gamm1$edf)) >= 1.99,
                                         ifelse(as.numeric(as.character(summary.gamm1$delta.AIC)) >= 2.0, 
                                                "yes", 
                                                "no"),
                                         "no"),
                                  "no")

summary.lmac$best.model = ifelse(as.numeric(as.character(summary.gamm1$p.ac))<=0.05,
                                 ifelse(as.numeric(as.character(summary.gamm1$edf))>=1.99,
                                        ifelse(as.numeric(as.character(summary.gamm1$delta.AIC))>=2.0,"no","yes"),"yes"),"no")

summary.gam1$best.model = ifelse(as.numeric(as.character(summary.gamm1$p.ac))>0.05,
                                 ifelse(as.numeric(as.character(summary.gam1$edf))>1.99,
                                        ifelse(as.numeric(as.character(summary.gam1$delta.GCV))<0,
                                               ifelse(as.numeric(as.character(summary.gam1$delta.AIC))<=-2.0,"yes","no"),"no"),"no"),"no")

summary.linear$best.model = ifelse(as.numeric(as.character(summary.gamm1$p.ac))>0.05,
                                   ifelse(as.numeric(as.character(summary.gam1$edf))<1.99,"yes",
                                          ifelse(as.numeric(as.character(summary.gam1$delta.GCV))>0,"yes",
                                                 ifelse(as.numeric(as.character(summary.gam1$delta.AIC))<=-2.0,"no","yes"))),"no")


allSummary <- rbind(summary.gamm1,
                    summary.lmac,
                    summary.gam1,
                    summary.linear)

subset(allSummary, best.model == "yes")

pick.gam <- which(allSummary$best.model=="yes" & allSummary$MODEL=="GAM")

####################################################################################################################
###CODE FOR MAKING PLOTs FOR RELATIONSHIPS THAT SHOWED GAM WAS THE BEST MODEL (using "pick.gam" from line 353 above)
####################################################################################################################

###############################################################################################################
## Changes from the GAMM bootstrapping and plotting above include:
## 1)Using gam instead of gamm
## 2)Remove gaussian bootstrap process, but use the resids from myresponse~gam1$fitted (this is different from Scott's original code). Using the resids gives CI around the gam fit that look more believable.
###############################################################################################################

gam1 <- gam(myresponse ~ s(mydriver, bs = "tp", k = ks), method = "GCV.Cp", se = T)

ind.fit <- list(mydriver = seq(from = min(mydriver), to = max(mydriver), length.out = sp.len))
pred <- predict.gam(gam1, newdata = ind.fit, type = "response", se.fit = T)

gam.data <- data.frame(pred$fit,
                       ind.fit$mydriver)
colnames(gam.data) <- c("response","driver")
gam.check(gam1)

# Bootstrap and subsequent code in a nutshell...
# 1) take a random sample of the data 1000 times
# 2) take GAM of each bootstrap replicate (br)
# 3) Take the 1st and 2nd derivative of each br
# 4) When 95% of the 1st or 2nd derivative replicates are greater or less than zero, we have a significant
#    trend (1st derivative), or threshold (2nd derivative).
#
# Note on maximum entropy bootstrap:
# Bootstrapping is a nonparametric alternative for assessing the 
# uncertainty of a linear trend in the presence
# of autocorrelation. Maximum entropy bootstrap allows
# us to construct a set of replicates of the original series
# to be used for inference while retaining the temporal
# dependence structure of the original series in the resampling process (Kirk & Stumpf 2009)


# respboot <- meboot(myresponse, reps = nb, trim = 0.1)$ens 
# drivboot <- meboot(mydriver, reps = nb, trim = 0.1)$ens
resids <- myresponse - gam1$fitted
plot(myresponse, resids)
###This is an attempt to reproduce the gaussian process regression bootstrap for the GAM...NOT COMPLETED YET
#dist <- diag(length(resids))
#for (i in 1:nrow(dist)) for (j in 1:ncol(dist)) dist[i,j] = abs(j-i)
#resids.use <- rmvnorm(nb,rep(mean(resids),length(resids)),sigma=var(resids)*dist,method="chol")
#resids.use <- rnorm(nb,rep(mean(resids),length(resids)),sd=sd(resids)*dist)

###The following is the code used above in the GAMM section
# sample multivariate normal for use in gaussian process regression bootstrap
#the gaussian bootstrap generates multivariate normal from the mean and a variance-covariance matrix of the residuals, where the covariance of the residuals is the residual variance multiplied by a correlation matrix based on the autocorrelation from the call to lme

#phi = coef(gam1$lme$modelStruct$corStruct,unconstrained=FALSE) #these coefficients come from a model using the "gamm" call...not a "gam" call
#dist <- matrix(0,nrow=length(resids),ncol=length(resids))
#for (i in 1:nrow(dist)) for (j in 1:ncol(dist)) dist[i,j] = abs(j-i)
#dist <- phi^dist
#library(mvtnorm)
#resids.use <- rmvnorm(nb,rep(mean(resids),length(resids)),sigma=var(resids)*dist,method="chol")

respboot <- matrix(data = NA, ncol = nb, nrow = length(mydriver))
#Use the following code for residual bootstrapping 
for (i in 1:nb) {
  #use next line if just doing residual bootstrap
  respboot[,i] <- gam1$fitted + sample(resids, length(resids), replace =TRUE)
  #use next line if doing gaussian process residual bootstrap
  #respboot[,i] <- gam1$fitted + resids.use[i,]
  drivboot[,i] <- mydriver
}


gam.mat <- matrix(ncol = nb, nrow = sp.len)
for(i in 1:nb) {
  resp.i <- respboot[,i]
  driv.i <- drivboot[,i]
  gam.i <- gam(resp.i ~ s(driv.i, bs = "tp", k = ks), method = "GCV.Cp", se = T)
  newdata.i <- list(driv.i = seq(from = min(driv.i), 
                                 to = max(driv.i),
                                 length.out = sp.len))
  gam.mat[,i] <- predict.gam(gam.i, newdata = newdata.i, type = "response", se.fit = T)$fit
}

#######################
# CI of GAM bootstrap #
#######################
ci <- matrix(nrow = 2, ncol = sp.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci) <- c("lower","upper")
for(j in 1:sp.len) {
  IC <- quantile(gam.mat[j,], c(0.025, 0.975), na.rm = TRUE)
  ci[,j] <- rbind(IC[1], IC[2])
}

#############################
# 1st Derivative Estimation #
#############################
dif1.line <- diff(gam.data$response, difference = 1) # Actual 1st deriv estimate from original smoother
deriv.matrix.1 <- matrix(nrow = sp.len -dif1, ncol = nb) ## create a matrix to for the 1st deriv estimates
for(k in 1:nb) {
  derivi <- gam.mat[,k]
  deriv.response <- diff(derivi, difference = dif1)
  driver.len <- length(mydriver) - dif1
  deriv.object <- cbind(deriv.response, driver.len)
  deriv.matrix.1[,k] <- deriv.object[,1]
}


# CI of 1st derivative 
dif1.len <- sp.len - dif1
ci1 <- matrix(nrow = 2, ncol = dif1.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci1) <- c("lower","upper")
#
for(l in 1:dif1.len) {
  IC<-quantile(deriv.matrix.1[l,], c(0.025, 0.975))
  ci1[,l]<-rbind(IC[1], IC[2])
}
#   
#############################
# 2nd Derivative Estimation #
#############################
dif2.line<-diff(gam.data$response, difference=2) # Actual 2nd deriv estimate from original smoother
deriv.matrix.2<-matrix(nrow=sp.len-dif2, ncol=nb) ## create a matrix to for the 2nd deriv estimates
for(m in 1:nb) {
  derivi<-gam.mat[,m]
  deriv.response<-diff(derivi, difference=dif2)
  driver.len<-length(mydriver)-dif2
  deriv.object <- cbind(deriv.response,driver.len)
  deriv.matrix.2[,m]<-deriv.object[,1]
}

# CI of 2nd derivative 
dif2.len<-sp.len-dif2
ci2<-matrix(nrow= 2, ncol= dif2.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci2)<-c("lower","upper")
#
for(n in 1:dif2.len) {
  IC<-quantile(deriv.matrix.2[n,], c(0.025, 0.975))
  ci2[,n]<-rbind(IC[1], IC[2])
}
#
# CI of response
lower <- min(ci[1,  ])
upper <- max(ci[2,  ])
lower1 <- min(ci1[1,  ])
upper1 <- max(ci1[2,  ])
lower2 <- min(ci2[1,  ])
upper2 <- max(ci2[2,  ])
#
CI.table<-matrix(nrow=2, ncol=3, data=c(upper,
                                        lower,
                                        upper1,
                                        lower1,
                                        upper2,
                                        lower2),
                 dimnames=list(c("UPPER", "LOWER"),
                               c("GAM","FIRST","SECOND")))          
CI.mat <- data.frame(CI.table)  #


###################
# Make some plots #
###################
#
changepts.lower1 <- seq(1, length(dif1.line)) [ci1["lower",  ] > 0]
changepts.upper1 <- seq(1, length(dif1.line)) [ci1["upper",  ] < 0]
changepts.lower2 <- seq(1, length(dif2.line)) [ci2["lower",  ] > 0]
changepts.upper2 <- seq(1, length(dif2.line)) [ci2["upper",  ] < 0]
#
driver1 <- gam.data$driver[-1]
driver2 <- gam.data$driver[-c(1,2)]
#

####################################################################################
## PLOTS function with 2nd derivative arrow in one-panel plot##
####################################################################################

# Will make a plot in your "figures directory"

#set your desktop 
# setwd("~/Desktop")
png("driver_response.png",width = 5, height = 3, units = "in", res = 600)

#########################################
# Figure: GAM pressure-response plot #
#########################################
par(mfg= c(1,1), 
    mar=c(2.75,2.75,0.5,0.25), 
    mgp=c(0,0.25,0), 
    cex= 0.75, 
    tck=-0.02, 
    family = "sans")

plot(myresponse ~ mydriver,
     #ylim = c(min(CI.mat$GAM[2],myresponse),max(CI.mat$GAM[1],myresponse)),
     ylim = c(ifelse(CI.mat$GAM[2] < min(myresponse), CI.mat$GAM[2], min(myresponse)),
              ifelse(CI.mat$GAM[1] > max(myresponse), CI.mat$GAM[1], max(myresponse))),
     type = "n",
     ylab="",
     xlab="",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 6)
axis(1, dri.ax,cex=1)
mtext("mydriver", side = 1, line = 1.5, cex=1)
#
gam.ax <- pretty(c(CI.mat$GAM[2], CI.mat$GAM[1]), 3)
axis(2, gam.ax,cex=1)
mtext("myresponse", side= 2,line=1.5, cex=1)
#
# Shading 95% CIs
polygon(c(gam.data$driver, rev(gam.data$driver)), 
        c(ci[2,  ], rev(ci[1,  ])),
        col = gray(0.95), 
        border = gray(.80))
#
# GAM Line
lines(gam.data$response ~ gam.data$driver, lty =2,  lwd=2)


# PLOT DERIVATIVES
points(gam.data$response[ci2["lower",  ] > 0] ~ gam.data$driver[ci2["lower",  ] > 0], 
       col = gray(.35), pch = 16)
points(gam.data$response[ci2["upper",  ] < 0] ~ gam.data$driver[ci2["upper",  ] < 0], 
       col = gray(.35), pch = 16)

# Add original data points
points(myresponse ~ mydriver, pch = 16, cex = 1, lwd = 1)


#Add dotted line at most probable threshold level of pressure
threshold.df <- data.frame(response = gam.data$response[-c(1,2)], driver2=driver2, dif2=dif2.line, 
                           ci2.upper=ci2["upper",], ci2.lower=ci2["lower",])
threshold.df.upper <- threshold.df[threshold.df$ci2.upper < 0, ]
threshold.upper.mlt <- data.frame(driver2 = threshold.df.upper$driver2[which.min(threshold.df.upper$dif2)],
                                  response = threshold.df.upper$response[which.min(threshold.df.upper$dif2)])
p=length(row.names(threshold.df.upper))  
ylim = c(min(CI.mat$GAM[2],myresponse),max(CI.mat$GAM[1],myresponse))
ifelse(ylim[1] > 0,
       b <- ylim[1]-((ylim[2]-ylim[1])*0.04),
       b <- ylim[1]*1.08)
threshold.df.lower <- threshold.df[threshold.df$ci2.lower > 0, ]
threshold.lower.mlt <- data.frame(driver2 = threshold.df.lower$driver2[which.max(threshold.df.lower$dif2)],
                                  response = threshold.df.lower$response[which.max(threshold.df.lower$dif2)])
threshold.lower.start <- threshold.df.lower[1,2]
threshold.lower.end <- threshold.df.lower[dim(threshold.df.lower)[[1]],2]
q = length(row.names(threshold.df.lower))

if(p!=0)segments(threshold.upper.mlt[[1]],b,threshold.upper.mlt[[1]],threshold.upper.mlt[[2]], lty=3, col = "red", lwd=2)
if(q!=0)segments(threshold.lower.mlt[[1]],b,threshold.lower.mlt[[1]],threshold.lower.mlt[[2]], lty=3, col = "red", lwd=2)

#Add arrows pointing to threshold level of pressure
arr.length <- 0.01*(gam.ax[4] - gam.ax[1])         #make arrow 10% the length of the y-axis
arrhead.length <- 0.1    # make arrowhead 
if(p!=0)arrows(x0 = threshold.upper.mlt[[1]], x1 = threshold.upper.mlt[[1]], y0 = b + arr.length,
               y1 = b,length=arrhead.length, col = "red", lwd=2)
if(q!=0)arrows(x0 = threshold.lower.mlt[[1]], x1 = threshold.lower.mlt[[1]], y0 = b + arr.length,
               y1 = b,length=arrhead.length, col = "red", lwd=2)

box()
dev.off()
##

####################################################################################
## PLOTS function, 1st derivative and 2nd derivative in multi-panel plot##
####################################################################################
# Will make a plot in your "figures directory"
# der.plot <- paste0(results.dir, paste0(dataset,".",resp.name, "-", "mydriver", ".3.panel.gam.png"))
# png(file = der.plot, width = 90, height = 130, units = "mm", res = 600)
# layout(matrix(c(1:3), 3,1,  byrow = T),
#        widths = lcm(8.45), heights = lcm(c(4.95,3.45, 4.45)),
#        respect = F)
#

png("driver_derivatives.png",width = 3, height = 8, units = "in", res = 600)


#########################################
# Figure A): GAM pressure-response plot #
#########################################
par(mfrow=c(3,1),
    mfg= c(1,1), 
    mar=c(0,2.15,0.5,0.25), 
    mgp=c(0,0.25,0), 
    cex= 0.75, 
    tck=-0.02, 
    family = "sans")

plot(myresponse ~ mydriver,
     #ylim = c(min(CI.mat$GAM[2],myresponse),max(CI.mat$GAM[1],myresponse)),
     ylim = c(ifelse(CI.mat$GAM[2] < min(myresponse), CI.mat$GAM[2], min(myresponse)),
              ifelse(CI.mat$GAM[1] > max(myresponse), CI.mat$GAM[1], max(myresponse))),
     type = "n",
     ylab="",
     xlab="",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 6)
axis(1, dri.ax, labels=F)
#mtext(dri.name, side = 1, line = 1.25, cex = 0.75)
#
gam.ax <- pretty(c(CI.mat$GAM[2], CI.mat$GAM[1]), 3)
axis(2, gam.ax)
mtext(resp.name, side= 2,line=1.5, cex= 0.75)
#
# Shading 95% CIs
polygon(c(gam.data$driver, rev(gam.data$driver)), 
        c(ci[2,  ], rev(ci[1,  ])),
        col = gray(0.95), 
        border = gray(.80))
#
# GAM Line
lines(gam.data$response ~ gam.data$driver, lty =2,  lwd=2)


# PLOT DERIVATIVES
points(gam.data$response[ci2["lower",  ] > 0] ~ gam.data$driver[ci2["lower",  ] > 0], 
       col = gray(.35), pch = 16)
points(gam.data$response[ci2["upper",  ] < 0] ~ gam.data$driver[ci2["upper",  ] < 0], 
       col = gray(.35), pch = 16)

# Add original data points
points(myresponse~mydriver, pch = 16, cex=1, lwd=1)
#text(myresponse~mydriver, labels = as.character(ts.length), cex=.5, offset = 0.6)  #adds points with year text
#

#Add dotted line at most probable threshold level of pressure
threshold.df <- data.frame(response = gam.data$response[-c(1,2)], driver2=driver2, dif2=dif2.line, 
                           ci2.upper=ci2["upper",], ci2.lower=ci2["lower",])
threshold.df.upper <- threshold.df[threshold.df$ci2.upper < 0, ]
threshold.upper.mlt <- data.frame(driver2 = threshold.df.upper$driver2[which.min(threshold.df.upper$dif2)],
                                  response = threshold.df.upper$response[which.min(threshold.df.upper$dif2)])
p=length(row.names(threshold.df.upper))  
ylim = c(min(CI.mat$GAM[2],myresponse),max(CI.mat$GAM[1],myresponse))
ifelse(ylim[1] > 0,
       b <- ylim[1]-((ylim[2]-ylim[1])*0.04),
       b <- ylim[1]*1.08)
threshold.df.lower <- threshold.df[threshold.df$ci2.lower > 0, ]
threshold.lower.mlt <- data.frame(driver2 = threshold.df.lower$driver2[which.max(threshold.df.lower$dif2)],
                                  response = threshold.df.lower$response[which.max(threshold.df.lower$dif2)])
threshold.lower.start <- threshold.df.lower[1,2]
threshold.lower.end <- threshold.df.lower[dim(threshold.df.lower)[[1]],2]
q=length(row.names(threshold.df.lower))  
if(p!=0)segments(threshold.upper.mlt[[1]],b,threshold.upper.mlt[[1]],threshold.upper.mlt[[2]], lty=3, col = "red", lwd=2)
if(q!=0)segments(threshold.lower.mlt[[1]],b,threshold.lower.mlt[[1]],threshold.lower.mlt[[2]], lty=3, col = "red", lwd=2)

#Add arrows pointing to threshold level of pressure
arr.length <- 0.01*(gam.ax[4] - gam.ax[1])         #make arrow 10% the length of the y-axis
arrhead.length <- 0.1    #make arrowhead 
if(p!=0)arrows(x0 = threshold.upper.mlt[[1]], x1 = threshold.upper.mlt[[1]], y0 = b + arr.length,
               y1 = b,length=arrhead.length, col = "red", lwd=2)
if(q!=0)arrows(x0 = threshold.lower.mlt[[1]], x1 = threshold.lower.mlt[[1]], y0 = b + arr.length,
               y1 = b,length=arrhead.length, col = "red", lwd=2)

# # Add legend
# legend("topleft",
#        legend="a)",
#        bty = "n",
#        box.lwd=0, 
#        box.col="white")
# #
# box()

########################
# Figure B) s'(X) plot #
########################
par(mfg = c(2,1), 
    mar=c(0,2.15,1,0.25), 
    mgp=c(0,0.25,0), 
    cex= 0.75, 
    tck=-0.02, 
    family = "sans")
#
plot(dif1.line ~ driver1,
     ylim =  c(CI.mat$FIRST[2], CI.mat$FIRST[1]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 6)
axis(1, dri.ax, tck=-0.02,labels=F)
#
d1.ax <- pretty(c(-.85,.85)*max(abs(CI.mat$FIRST)), 6)
axis(2, d1.ax)
mtext("s'(X)",side= 2,line=1.5, cex= .75)

#This polygon draws a line on the x-axis where the upper CI is < 0
polygon(c(driver1, rev(driver1)), 
        c(ci1[2,], rev(ci1[1,])),
        col = gray(0.95), border = gray(.80))

#This polygon draws a line on the x-axis where the upper CI is < 0
polygon(c(driver1[changepts.upper1], rev(driver1[changepts.upper1])),
        c(ci1[2,][changepts.upper1], rev(as.matrix(rep(max(ci1[2,][changepts.upper1]),
                                                       length(changepts.upper1)))[,1])),col="black", border="black",lwd=2)

#This polygon draws a line on the x-axis where the lower CI is > 0
polygon(c(driver1[changepts.lower1], rev(driver1[changepts.lower1])),
        c(ci1[1,][changepts.lower1], rev(as.matrix(rep(min(ci1[1,][changepts.lower1]),
                                                       length(changepts.lower1)))[,1])),col="black", border="black",lwd=2)
#
abline(h=0, lwd=.5)
lines(dif1.line ~ driver1, lwd=1.5)
#
# Arrow stuff should be updated ad hoc
# arrows(x0=gam.data$driver[changepts.lower1][4], y0= min(ci1[1,])*.75,
#        x1=gam.data$driver[changepts.lower1][4], y1= min(ci1[1,])*.15, 
#        col="black", 
#        length=.10,
#        angle=25,
#        lwd=1)
# 
# legend("topleft",
#        legend="b)",
#        bty = "n",
#        box.lwd=0, 
#        box.col="white")
# #
# box()
########################
# Figure C) s"(X) plot #
########################
########################
# Figure C) s"(X) plot #
########################

par(mfg = c(3,1), 
    mar=c(2.15,2.15,1,0.25), 
    mgp=c(0,0.25,0), 
    cex= 0.75, 
    tck=-0.02, 
    family = "sans")
#
plot(dif2.line ~ driver2,
     ylim =  c(CI.mat$SECOND[2], CI.mat$SECOND[1]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 6)
axis(1, tck=-0.02,dri.ax)
mtext("mydriver", side = 1, line = 1.25, cex=0.75)
#
d2.ax <- pretty(c(-.85,.85) * max(abs(CI.mat$SECOND)), 3)
axis(2, d2.ax)
mtext('s"(X)', side= 2,line=1.5, cex= .75)
#
#This polygon draws the CI around the 2nd derivative line
polygon(c(driver2, rev(driver2)), 
        c(ci2[2,], rev(ci2[1,])),
        col = gray(0.95), border = gray(.80))

#This polygon draws a line on the x-axis where the upper CI is < 0
polygon(c(driver2[changepts.upper2], rev(driver2[changepts.upper2])),
        c(ci2[2,][changepts.upper2], rev(as.matrix(rep(max(ci2[2,][changepts.upper2]), 
                                                       length(changepts.upper2)))[,1])),col="black", border="black", lwd=2)

#This polygon draws a line on the x-axis where the lower CI is > 0
polygon(c(driver2[changepts.lower2], rev(driver2[changepts.lower2])), 
        c(ci2[1,][changepts.lower2], rev(as.matrix(rep(min(ci2[1,][changepts.lower2]), 
                                                       length(changepts.lower2)))[,1])),col="black", border="black",lwd=2)

#
abline(h=0, lwd=.5)
lines(dif2.line ~ driver2, lwd=1.5)
# 
# legend("topleft",
#        legend="c)",
#        bty = "n",
#        box.lwd = 0, 
#        box.col = "white")



dev.off()



###############################################################
#Magnitude of change calculations (used to measure effect size)
###############################################################
#Below, "lower" refers to thresholds in which the lower CI of the 2nd derivative is >0 and "upper" refers to thresholds in which the upper CI is <0; "mlt" refers to the "most likely threshold value where the 2nd derivative is at its greatest value (positive or negative) within the threshold range; "start" refers to the first driver value in the threshold range and "end" refers to the last driver value in the threshold range;"1" refers to the mean of the gam line to the left of threshold, "2" refers to the mean of the gam line to the right of threshold.

#1.Identify 3 points in threshold(s) range: start, maximum value of 2nd derivative within range, and end
start.lower <- ifelse(q!=0,threshold.df.lower[1,2],NA)
mlt.lower <- ifelse(q!=0,threshold.lower.mlt[[1]],NA) #most likely threshold (mlt)
end.lower <- ifelse(q!=0,threshold.df.lower[dim(threshold.df.lower)[[1]],2],NA)

start.upper <- ifelse(p!=0,threshold.df.upper[1,2],NA)
mlt.upper <- ifelse(p!=0,threshold.upper.mlt[[1]],NA) #most likely threshold (mlt)
end.upper <- ifelse(p!=0,threshold.df.upper[dim(threshold.df.upper)[[1]],2],NA)

#2.Calculate response means and medians on both sides of threshold points using gam predicted values (i.e. gam line)
start.mean.1.lower <- mean(gam.data$response[which(gam.data$driver<start.lower)])
start.mean.2.lower <- mean(gam.data$response[which(gam.data$driver>start.lower)])
mlt.mean.1.lower <- mean(gam.data$response[which(gam.data$driver<mlt.lower)])
mlt.mean.2.lower <- mean(gam.data$response[which(gam.data$driver>mlt.lower)])
end.mean.1.lower <- mean(gam.data$response[which(gam.data$driver<end.lower)])
end.mean.2.lower <- mean(gam.data$response[which(gam.data$driver>end.lower)])

start.mean.1.upper <- mean(gam.data$response[which(gam.data$driver<start.upper)])
start.mean.2.upper <- mean(gam.data$response[which(gam.data$driver>start.upper)])
mlt.mean.1.upper <- mean(gam.data$response[which(gam.data$driver<mlt.upper)])
mlt.mean.2.upper <- mean(gam.data$response[which(gam.data$driver>mlt.upper)])
end.mean.1.upper <- mean(gam.data$response[which(gam.data$driver<end.upper)])
end.mean.2.upper <- mean(gam.data$response[which(gam.data$driver>end.upper)])

start.median.1.lower <- median(gam.data$response[which(gam.data$driver<start.lower)])
start.median.2.lower <- median(gam.data$response[which(gam.data$driver>start.lower)])
mlt.median.1.lower <- median(gam.data$response[which(gam.data$driver<mlt.lower)])
mlt.median.2.lower <- median(gam.data$response[which(gam.data$driver>mlt.lower)])
end.median.1.lower <- median(gam.data$response[which(gam.data$driver<end.lower)])
end.median.2.lower <- median(gam.data$response[which(gam.data$driver>end.lower)])

start.median.1.upper <- median(gam.data$response[which(gam.data$driver<start.upper)])
start.median.2.upper <- median(gam.data$response[which(gam.data$driver>start.upper)])
mlt.median.1.upper <- median(gam.data$response[which(gam.data$driver<mlt.upper)])
mlt.median.2.upper <- median(gam.data$response[which(gam.data$driver>mlt.upper)])
end.median.1.upper <- median(gam.data$response[which(gam.data$driver<end.upper)])
end.median.2.upper <- median(gam.data$response[which(gam.data$driver>end.upper)])

#3.Calculate proportional change on either side of threshold points
start.diff.mean.lower <- start.mean.2.lower - start.mean.1.lower
mlt.diff.mean.lower <- mlt.mean.2.lower - mlt.mean.1.lower
end.diff.mean.lower <- end.mean.2.lower - end.mean.1.lower
#
start.mag.mean.lower <- ifelse(q!=0,start.diff.mean.lower/abs(start.mean.1.lower),NA)
mlt.mag.mean.lower <- ifelse(q!=0,mlt.diff.mean.lower/abs(mlt.mean.1.lower),NA)
end.mag.mean.lower <- ifelse(q!=0,end.diff.mean.lower/abs(end.mean.1.lower),NA)
####
start.diff.mean.upper <- start.mean.2.upper - start.mean.1.upper
mlt.diff.mean.upper <- mlt.mean.2.upper - mlt.mean.1.upper
end.diff.mean.upper <- end.mean.2.upper - end.mean.1.upper
#
start.mag.mean.upper <- ifelse(p!=0,start.diff.mean.upper/abs(start.mean.1.upper),NA)
mlt.mag.mean.upper <- ifelse(p!=0,mlt.diff.mean.upper/abs(mlt.mean.1.upper),NA)
end.mag.mean.upper <- ifelse(p!=0,end.diff.mean.upper/abs(end.mean.1.upper),NA)
#
start.diff.median.lower <- start.median.2.lower - start.median.1.lower
mlt.diff.median.lower <- mlt.median.2.lower - mlt.median.1.lower
end.diff.median.lower <- end.median.2.lower - end.median.1.lower
#
start.mag.median.lower <- ifelse(q!=0,start.diff.median.lower/abs(start.median.1.lower),NA)
mlt.mag.median.lower <- ifelse(q!=0,mlt.diff.median.lower/abs(mlt.median.1.lower),NA)
end.mag.median.lower <- ifelse(q!=0,end.diff.median.lower/abs(end.median.1.lower),NA)
####
start.diff.median.upper <- start.median.2.upper - start.median.1.upper
mlt.diff.median.upper <- mlt.median.2.upper - mlt.median.1.upper
end.diff.median.upper <- end.median.2.upper - end.median.1.upper
#
start.mag.median.upper <- ifelse(p!=0,start.diff.median.upper/abs(start.median.1.upper),NA)
mlt.mag.median.upper <- ifelse(p!=0,mlt.diff.median.upper/abs(mlt.median.1.upper),NA)
end.mag.median.upper <- ifelse(p!=0,end.diff.median.upper/abs(end.median.1.upper),NA)
#

resp.start.lower <- ifelse(q!=0,threshold.df.lower[[1,1]],NA)
resp.mlt.lower <- ifelse(q!=0,threshold.lower.mlt[[2]],NA)
resp.end.lower <- ifelse(q!=0,threshold.df.lower[[dim(threshold.df.lower)[1],1]],NA)
dri.start.lower <- ifelse(q!=0,threshold.df.lower[[1,2]],NA)
dri.end.lower <- ifelse(q!=0,threshold.df.lower[[dim(threshold.df.lower)[1],2]],NA)
resp.start.upper <- ifelse(p!=0,threshold.df.upper[[1,1]],NA)
resp.mlt.upper <- ifelse(p!=0,threshold.upper.mlt[[2]],NA)
resp.end.upper <- ifelse(p!=0,threshold.df.upper[[dim(threshold.df.upper)[1],1]],NA)
dri.start.upper <- ifelse(p!=0,threshold.df.upper[[1,2]],NA)
dri.end.upper <- ifelse(p!=0,threshold.df.upper[[dim(threshold.df.upper)[1],2]], NA)

magnitude <- as.data.frame(cbind("myresponse",
                                 "mydriver",
                                 resp.start.lower, #response value at beginning of threshold
                                 resp.mlt.lower, #response value at most likely threshold value
                                 resp.end.lower, #response value at end of threshold
                                 start.lower, #driver value at beginning of threshold
                                 mlt.lower, #driver value where 2nd derivative is furthest away from 0
                                 end.lower, #driver value at end of threshold
                                 start.mag.mean.lower, #magnitude of response across the threshold at start of threshold range
                                 mlt.mag.mean.lower, #magnitude of response across the threshold at most likely threshold value of threshold range
                                 end.mag.mean.lower, #magnitude of response across the threshold at end of threshold range
                                 start.mag.median.lower, #magnitude of response across the threshold at start of threshold range
                                 mlt.mag.median.lower, #magnitude of response across the threshold at most likely threshold value of threshold range
                                 end.mag.median.lower, #magnitude of response across the threshold at end of threshold range
                                 resp.start.upper, #response value at beginning of threshold
                                 resp.mlt.upper, #response value at most likely threshold value
                                 resp.end.upper, #response value at end of threshold
                                 start.upper, #driver value at beginning of threshold
                                 mlt.upper, #driver value where 2nd derivative is furthest away from 0
                                 end.upper, #driver value at end of threshold
                                 start.mag.mean.upper, #magnitude of response across the threshold at start of threshold range
                                 mlt.mag.mean.upper, #magnitude of response across the threshold at most likely threshold value of threshold range
                                 end.mag.mean.upper, #magnitude of response across the threshold at end of threshold range
                                 start.mag.median.upper, #magnitude of response across the threshold at start of threshold range
                                 mlt.mag.median.upper, #magnitude of response across the threshold at most likely threshold value of threshold range
                                 end.mag.median.upper #magnitude of response across the threshold at end of threshold range
                                 
))

colnames(magnitude)<-c("Response","Driver","Response.start.lower","Response.mlt.lower","Response.end.lower","Threshold.start.lower","Most.likely.threshold.lower","Threshold.end.lower","Mean magnitude.start.lower","Mean magnitude.mlt.lower","Mean magnitude.end.lower","Median magnitude.start.lower","Median magnitude.mlt.lower","Median magnitude.end.lower","Response.start.upper","Response.mlt.upper","Response.end.upper","Threshold.start.upper","Most.likely.threshold.upper","Threshold.end.upper","Mean magnitude.start.upper","Mean magnitude.mlt.upper","Mean magnitude.end.upper","Median magnitude.start.upper","Median magnitude.mlt.upper","Median magnitude.end.upper")

allMagnitude <- rbind(magnitude)

write.csv(allMagnitude,"magnitude of threshold.gam.csv",row.names=F)


} #end of test loop