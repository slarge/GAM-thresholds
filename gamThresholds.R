rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: Thresholds from a GAM using maximum entropy bootstrapping
# AUTHOR: Scott Large 2013 (revised 2015)
# REVIEWED BY:
# VERSION: 0.1
#

######################
# CHANGES/ ADDITIONS #
######################

# Need to add: 

# Done:

############
# PACKAGES #
############
library(meboot)
library(mgcv)
library(reshape2)
#
##################
# Set up folders #
##################
data.dir <- "data/"
figure.dir <- "output/"
#
######################
# Load any functions #
######################
#
#
#############
# Load data #
#############
set.seed(627)
cc <- read.csv("data/CCIEA-RPW.csv")
#
# Subset area... initally for coastwide
ccALL <- cc[cc$Coastwide == 1,]
ccALL$year <- as.numeric(ccALL$year)
#
ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
#
# Wide df with columns as variables
dat.full <- dcast(ccALL, year ~ timeseries, value.var = "value")
#
resp.name <- "GF_MTL_coastwide" 
dri.name <- "Coastal_engineering"
#
dat.cut <- dat.full[colnames(dat.full) %in% c("year", resp.name, dri.name)]
dat.cut <- dat.cut[apply(dat.cut, 1, function(x)!any(is.na(x))),]
#
myresponse <- dat.cut[[resp.name]]
mydriver <- dat.cut[[dri.name]]
ts.length <- dat.cut$year
#
#####################
# Fit the GAM model #
#####################
sp.len <- 200 # Spline length
nb <- 1000  # Number of bootstrap replicates
gam.mat <- matrix(nrow=sp.len, ncol=nb) ## create a matrix to be filled with bootstrapped splines (200 x 1000)
dif1 <- 1 # First derivative
dif2 <- 2 # Second derivative
#
# GAM #
gam1 <- gam(myresponse ~ s(mydriver, bs = "ts"), method = "GCV.Cp", se = T)
# gam1 uses thin-plate splines and prevents overfitting by using cross validation, run 
# ?gam() for more info on model. 
# summary(gam1) will show the summary statistics for the model
# gam.check(gam1)
# plot(gam1) will show a plot of the model

# GLM #
gam2 <- gam(myresponse ~ mydriver, method = "GCV.Cp", se = T)
# gam2 removes the smoothing function s(), so it is basically a linear model. Check out my paper to see the conditions
# on when to pick a GAM over a GLM... although for your data I imagine that GAM will be best.
# gam.check(gam2)

# Pull out relevant model outputs:
summary.gam1 <- as.data.frame(cbind(summary(gam1)$dev.expl,       # deviance explained
                                    summary(gam1)$edf,            # estimated degrees of freedom
                                    summary(gam1)$sp.criterion,   # GCV score
                                    summary(gam1)$s.pv))          # estimated p-value
colnames(summary.gam1)<-c("dev.expl","edf", "GCV", "p_value")

ind.fit <- list(mydriver = seq(from = min(mydriver), to = max(mydriver), length.out = sp.len))
pred <- predict.gam(gam1, newdata = ind.fit, type = "response", se.fit = T)
#
gam.data <- data.frame(pred$fit,
                       ind.fit$mydriver)
colnames(gam.data) <- c("response","driver")

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
# dependence structure of the original series in the resampling process
#
respboot <- meboot(myresponse, reps = nb, trim = 0.1)$ens 
drivboot <- meboot(mydriver, reps = nb, trim = 0.1)$ens
#
gam.mat <- matrix(ncol = nb, nrow = sp.len)
for(i in 1:nb) {
  resp.i <- respboot[,i]
  driv.i <- drivboot[,i]
  gam.i <- gam(resp.i ~ s(driv.i, bs = "ts"), method = "GCV.Cp", se = T)
  newdata.i <- list(driv.i = seq(from = min(driv.i), 
                                 to = max(driv.i),
                                 length.out = sp.len))
  gam.mat[,i] <- predict.gam(gam.i, newdata = newdata.i, type = "response", se.fit = T)$fit
}
#
#######################
# CI of GAM bootstrap #
#######################
ci <- matrix(nrow= 2, ncol= sp.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci)<-c("lower","upper")
for(j in 1:sp.len) {
  IC <- quantile(gam.mat[j,], c(0.025, 0.975))
  ci[,j]<-rbind(IC[1], IC[2])
}

#############################
# 1st Derivative Estimation #
#############################
dif1.line<-diff(gam.data$response, difference=1) # Actual 1st deriv estimate from original smoother
deriv.matrix.1<-matrix(nrow=sp.len-dif1, ncol=nb) ## create a matrix to for the 1st deriv estimates
for(k in 1:nb) {
  derivi<-gam.mat[,k]
  deriv.response<-diff(derivi, difference=dif1)
  driver.len<-length(mydriver)-dif1
  deriv.object <- cbind(deriv.response,driver.len)
  deriv.matrix.1[,k]<-deriv.object[,1]
}
  
   
# CI of 1st derivative 
dif1.len <- sp.len-dif1
ci1 <- matrix(nrow= 2, ncol= dif1.len) ## create a matrix to be filled with bootstrapped CI
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
CI.mat <- data.frame(CI.table)
#
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
###########
## PLOTS ##
###########
# Will make a plot in your "figures directory"
der.plot <- paste0(figure.dir,"CCE_CE-MTL_v01.png")
png(file = der.plot, width = 170, height = 130, units = "mm", res = 600)
layout(matrix(c(1:3), 3,1,  byrow = T),
       widths = lcm(8.45), heights = lcm(c(4.95,3.95, 3.95)),
       respect = F)
#
####################################
# Figure A) pressure-response plot #
####################################
par(mfg= c(1,1), 
    mar=c(0,2.15,0,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")

plot(myresponse ~ mydriver,
     ylim = c(CI.mat$GAM[2], CI.mat$GAM[1]),
     type = "n",
     ylab="",
     xlab="",
     axes = F)
#
axis(1, tck = .015, labels= F)
#
gam.ax <- pretty(c(CI.mat$GAM[2], CI.mat$GAM[1]), 3)
axis(2, gam.ax, cex=.75)
mtext(resp.name, side= 2,line=1.25, cex= .75)
#
# Shading 95% CIs
polygon(c(gam.data$driver, rev(gam.data$driver)), 
        c(ci[2,  ], rev(ci[1,  ])),
        col = gray(0.95), 
        border = gray(.80))
#
# GAM Line
lines(gam.data$response ~ gam.data$driver, lty =2,  lwd=2)
#
# PLOT DERIVATIVES
points(gam.data$response[ci1["lower",  ] > 0] ~ gam.data$driver[ci1["lower",  ] > 0], 
      col = gray(.25), pch = 16)
points(gam.data$response[ci1["upper",  ] < 0] ~ gam.data$driver[ci1["upper",  ] < 0], 
      col = gray(.25), pch = 16)
points(gam.data$response[ci2["lower",  ] > 0] ~ gam.data$driver[ci2["lower",  ] > 0], 
      col = gray(.35), pch = 16)
points(gam.data$response[ci2["upper",  ] < 0] ~ gam.data$driver[ci2["upper",  ] < 0], 
      col = gray(.35), pch = 16)
#
# Add original data points
points(myresponse~mydriver, pch= 1, cex=.75, lwd=1)
#
# Add legend
legend("topleft",
       legend="a",
       bty = "n",
       box.lwd=0, 
       box.col="white")
#
# Box it up
box()
########################
# Figure B) s'(X) plot #
########################
par(mfg = c(2,1), 
    mar=c(0 ,2.15,0.15,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")
#
plot(dif1.line ~ driver1,
     ylim =  c(CI.mat$FIRST[2], CI.mat$FIRST[1]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
axis(1, tck = .015, labels= F)
#
d1.ax <- pretty(c(-.85,.85)*max(abs(CI.mat$FIRST)), 3)
axis(2, d1.ax, cex=.75)
mtext("s'(X)",side= 2,line=1.25, cex= .75)
#
polygon(c(driver1[changepts.upper1], rev(driver1[changepts.upper1])),
        c(ci1[2,][changepts.upper1], rev(as.matrix(rep(0, length(changepts.upper1)))[,1])),
        col="black", border=NA)
#
polygon(c(driver1[changepts.lower1], rev(driver1[changepts.lower1])),
        c(ci1[1,][changepts.lower1], rev(as.matrix(rep(0, length(changepts.lower1)))[,1])),
        col="black", border=NA)
#
polygon(c(driver1, rev(driver1)), 
        c(ci1[2,], rev(ci1[1,])),
        col = gray(0.95), border = gray(.80))
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
legend("topleft",
       legend="b",
       bty = "n",
       box.lwd=0, 
       box.col="white")
#
box()
########################
# Figure C) s"(X) plot #
########################
par(mfg = c(3,1), 
    mar=c(2.15,2.15,0.15,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")
#
plot(dif2.line ~ driver2,
     ylim =  c(CI.mat$SECOND[2], CI.mat$SECOND[1]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 3)
axis(1, dri.ax, cex = .75)
mtext(dri.name, side = 1, line = 1.25, cex = .75)
#
d2.ax <- pretty(c(-.85,.85) * max(abs(CI.mat$SECOND)), 3)
axis(2, d2.ax, cex=.75)
mtext('s"(X)', side= 2,line=1.25, cex= .75)
#
polygon(c(driver2[changepts.upper2], rev(driver2[changepts.upper2])),
        c(ci2[2,][changepts.upper2], rev(as.matrix(rep(0, length(changepts.upper2)))[,1])),
        col="black", border=NA)
#
polygon(c(driver2[changepts.lower1], rev(driver2[changepts.lower1])),
        c(ci2[1,][changepts.upper2], rev(as.matrix(rep(0, length(changepts.upper2)))[,1])),
        col="black", border=NA)
#
polygon(c(driver2, rev(driver2)), 
        c(ci2[2,], rev(ci2[1,])),
        col = gray(0.95), border = gray(.80))
#
abline(h=0, lwd=.5)
lines(dif2.line ~ driver2, lwd=1.5)
# 
legend("topleft",
       legend="c",
       bty = "n",
       box.lwd = 0, 
       box.col = "white")
#
box()
dev.off()
####