rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: Calculate the gradient of a multi-factor GAM
# AUTHOR: Scott Large 2013 (revised 2015)
# REVIEWED BY:
# VERSION: 0.1
#

######################
# CHANGES/ ADDITIONS #
######################
# Need to add: 
#
# Done:
#
############
# PACKAGES #
############
library(mgcv)
library(rgl)
library(meboot)
library(reshape2)
#
##################
# Set up folders #
##################
data.dir <- "data/"
figure.dir <- "output/"
source.dir <- "functions/"
#
######################
# Load any functions #
######################
source(paste0(source.dir, "bestGAM.R"))
source(paste0(source.dir, "myGrad.R"))
source(paste0(source.dir, "hullFunction.R"))
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
ind.name <- c("GF_spp_richness_coastwide", "GF-Simp_coastwide", "GF_MTL_coastwide", "Scav_ratio")
#
load(paste0(data.dir, "dfaTrends_v01.RDATA"))
#
driver <- melt(data.frame(ts.trends), id.vars = "YEAR",
               variable.name = "DRIVER",
               value.name = "VALUE")

indicator <- dat.full[colnames(dat.full) %in% c("year", ind.name)]
indicator <- melt(indicator, id.vars = "year",
                  variable.name = "INDICATOR", 
                  value.name = "VALUE")
colnames(indicator)[1] <- "YEAR"
#
#################
# Data Analysis #
#################

# Full model
ind.list <- unique(indicator$INDICATOR)
# 
# FM1 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3))        # LANDINGS
# FM2 <- formula(get(ei.name) ~ s(T1, k = fm.k3))              # T1
# FM3 <- formula(get(ei.name) ~ s(T2, k = fm.k3))              # T2
# # FM4 <- formula(get(ei.name) ~ s(T1, k = fm.k3) + s(T2, k = fm.k3))        # T1 + T2
# FM4 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3) + s(T1, k = fm.k3))  # LANDINGS + T1
# FM5 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3) + s(T2, k = fm.k3))  # LANDINGS + T2
# # FM7 <- formula(get(ei.name) ~ s(T1, k = fm.k3) + s(T2, k = fm.k3))        # LANDINGS + T1
# FM.list <- paste("FM", seq(1,5,1), sep = "")
# 

# Scaled value X original sd  + original mean 
# sv * sv.sd + sv.me

### INDIVIDUAL GAMS ####
test <- expand.grid(as.character(unique(indicator$INDICATOR)), 
                    c("T1", "T2", "T1"),
                    c("T2", "T3", "T3"))
colnames(test) <- c("INDICATOR", "TSA", "TSB")

test$TSA <- factor(test$TSA, levels = c("T1", "T2", "T3"))
test$TSB <- factor(test$TSB, levels = c("T1", "T2", "T3"))
#
gams.list <- list()
best.df <- data.frame()
mod.select <- data.frame()
#
for(ti in 1:nrow(test)) {
  bg.ti <-  test[ti,]
  bg.x <- driver[driver$DRIVER == bg.ti[,2],]
  bg.y <-  driver[driver$DRIVER == bg.ti[,3],]
  bg.z <-  indicator[indicator$INDICATOR == bg.ti[,1],]
  bg.name <- paste(bg.ti[,1], bg.ti[,2], bg.ti[,3], sep = ".")
  bg.year <- intersect(intersect(bg.x$YEAR, bg.y$YEAR), intersect(bg.x$YEAR, bg.z$YEAR))
  bg <- data.frame(x = scale(bg.x$VALUE[bg.x$YEAR %in% bg.year], scale = T, center = T),
                   y = scale(bg.y$VALUE[bg.y$YEAR %in% bg.year], scale = T, center = T),
                   z = scale(bg.z$VALUE[bg.z$YEAR %in% bg.year], scale = T, center = T),
                   YEAR = bg.year)
  bg <- bg[complete.cases(bg),]
  bg.year <- bg$YEAR
  bg$YEAR <- NULL
  if(length(bg$x) == length(bg$y)) {
    mat.len <- length(bg$x)
  } else stop("X and Y must be the same length")
  
  # Use GAM to model the relationship between these variables
  k.2 <- ceiling(length(bg.year)/2 - 2)
#   k.3 <- ceiling(length(bg.year)/4 - 2)
  
  F1 <- formula(z ~ s(x,  k = k.2))
  F2 <- formula(z ~ s(y,  k = k.2))
  F3 <- formula(z ~ s(x,  k = ceiling(k.2/2)) + s(y,  k = ceiling(k.2/2)))
  F4 <- formula(z ~ s(x,  k = k.2) + s(y,  k = k.2))
  F5 <- formula(z ~ x + s(y,  k = k.2))
  F6 <- formula(z ~ s(x,  k = k.2) + y)
  F7 <- formula(z ~ x + y)
  F8 <- formula(z ~ x)
  F9 <- formula(z ~ y)
  F.list <- paste("F", seq(1,9,1), sep = "")

  # Prints a list of all the model output
  gam.ti <- best.gam(formulaList = F.list, dataFrame = bg, gamMethod= "GCV.Cp")
  mod.select <- rbind(mod.select, cbind("MODEL.VARS" = bg.name, gam.ti))
    
  # Organize models based upon AICc
  best.aic <- gam.ti[order(gam.ti$AICc, decreasing = F),]
  
  # Make sure the gam() below matches best.gam(dataFrame, gamMethod)
  best.mod <- gam(get(best.aic[1,1]), method = "GCV.Cp", se = T, data = bg)
  best.df <- rbind(best.df, cbind("MODEL.VARS" = bg.name, best.aic[1,]))
  gams.list[[bg.name]] <- best.mod
  cat(paste("...", round(ti/nrow(test) * 100), "% /n", sep = ""))
}
save(best.df, file = "gamDFATable_v0008.RDATA")
write.csv(best.df, paste0(figure.dir, "dfaTable_v0008.csv"), row.names = F)
# 
     
# Bivariate GAM smoothers
fin.list <- best.df[best.df$X_EDF > 1 &
                      best.df$Y_EDF > 1 &
                      complete.cases(best.df),]
#
for(fl in 1:nrow(fin.list)) {
df.name <- as.character(fin.list[fl, 1])
gam1 <- gams.list[[df.name]]
#
df.parts <- unlist(strsplit(df.name, "[.]"))
df.x <- subset(driver, DRIVER == "LANDINGS")
df.y <- subset(driver, DRIVER == df.parts[2])
df.z <- subset(indicator, INDICATOR == df.parts[1])
#
df.year <- intersect(intersect(df.x$YEAR, df.y$YEAR), intersect(df.x$YEAR, df.z$YEAR))
df.ns <- data.frame(x =  df.x$VALUE[df.x$YEAR %in% df.year],
                    y =  df.y$VALUE[df.y$YEAR %in% df.year],
                    z =  df.z$VALUE[df.z$YEAR %in% df.year])
#
df <- data.frame(x = scale(df.x$VALUE[df.x$YEAR %in% df.year], scale = T, center = T),
                 y =  scale(df.y$VALUE[df.y$YEAR %in% df.year], scale = T, center = T),
                 z =  scale(df.z$VALUE[df.z$YEAR %in% df.year], scale = T, center = T))
#
# Create a convex hull surrounding data to aide in plotting
dia <- max(max(df$x) * .1, max(df$y) * .1)
hull <- hullFun(df.ns[,c(1:2)], diameter = dia)
#
if(length(df$x) == length(df$y)) {
  mat.len <- length(df$x)
} else stop("X and Y must be the same length")
#
######################
# Partial Derivative #
######################
# CONTROL PARAMETERS #
sp.len <- mat.len * mat.len # Spline length
# x.ran <- range(df$x)
# y.ran <- range(df$y)
x.ran <- range(hull$x)
y.ran <- range(hull$y)
eps <- 1e-7 # epsilon for finite difference
CV <- .95   # critical value for bootstrapped CI
set.seed(123)
nb <- 1000  # Bootstrap length
xseq <- seq(x.ran[1], to = x.ran[2], length.out = mat.len)
yseq <- seq(y.ran[1], to = y.ran[2], length.out = mat.len)
eps.list <- matrix(nrow = 9, ncol = 2, 
                   c(+eps, +eps, 0, 0, -eps, 0, -eps, +eps, -eps,
                     +eps, 0, +eps, 0, 0, -eps, -eps, -eps, +eps), 
                   byrow = F)
k.2 <- ceiling(length(df.year)/3 - 2)
k.3 <- ceiling(length(df.year)/4 - 2)
F.i <- summary(gam1)$formula # best.gam() formula
#
################### 
# Bootstrap setup #
###################
#######################
### Naive Bootstrap ###
#######################
# respboot <- matrix(replicate(nb, df$z[sample(mat.len, rep = TRUE)]), mat.len, nb)
# drivboot <- matrix(replicate(nb, df$x[sample(mat.len, rep = TRUE)]), mat.len, nb)
# enviboot <- matrix(replicate(nb, df$y[sample(mat.len, rep = TRUE)]), mat.len, nb)

################################################
### Maximum Entropy Bootstrapped Time Series ###
################################################
respboot <- meboot(df$z, reps = nb, trim = 0.1)$ens 
drivboot <- meboot(df$x, reps = nb, trim = 0.1)$ens
enviboot <- meboot(df$y, reps = nb, trim = 0.1)$ens

#############################
# Empty arrays for new data #
#############################
## Gradient
grad.Fx <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fy <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fxx <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fyy <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fxy <- array(NA, dim = c(mat.len, mat.len, nb))

grad.Fx.se <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fy.se <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fxx.se <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fyy.se <- array(NA, dim = c(mat.len, mat.len, nb))
grad.Fxy.se <- array(NA, dim = c(mat.len, mat.len, nb))

## X and Y range
dri.ran <- range(drivboot)
env.ran <- range(enviboot)

## X and Y sequence
dri.seq <- seq(dri.ran[1], to = dri.ran[2], length.out = mat.len)
env.seq <- seq(env.ran[1], to = env.ran[2], length.out = mat.len)

#####################
# BOOTSTRAP ROUTINE #
#####################
for(i in 1:nb) {
  # Isolate the data
  z <- respboot[,i]
  x <- drivboot[,i]
  y <- enviboot[,i]
  # Fit the GAM model
  gam.i <- gam(F.i,  method = "GCV.Cp", se = T)
  # calculate the finite differences for each component of the partial derivative
  tt <- my.grad(the.gam = gam.i, xs = dri.seq, ys = env.seq, eps.mat = eps.list)
  grad.Fx[,,i] <- matrix(tt$x$grad, mat.len, mat.len)  
  grad.Fy[,,i] <- matrix(tt$y$grad, mat.len, mat.len)
  grad.Fxx[,,i] <- matrix(tt$xx$grad, mat.len, mat.len)
  grad.Fyy[,,i] <- matrix(tt$yy$grad, mat.len, mat.len)
  grad.Fxy[,,i] <- matrix(tt$xy$grad, mat.len, mat.len)
  # calculate SE 
  grad.Fx.se[,,i] <- matrix(tt$x$se.grad, mat.len, mat.len)  
  grad.Fy.se[,,i] <- matrix(tt$y$se.grad, mat.len, mat.len)
  grad.Fxx.se[,,i] <- matrix(tt$xx$se.grad, mat.len, mat.len)  
  grad.Fyy.se[,,i] <- matrix(tt$yy$se.grad, mat.len, mat.len)
  grad.Fxy.se[,,i] <- matrix(tt$xy$se.grad, mat.len, mat.len)
  if(i %% 100 == 0) cat(paste("...", i/nb * 100, "% ", sep = ""))
}
# 

#######################
# CALCULATE QUANTILES #
#######################
# for each row and column of grad.Fx and grad.Fy take the quantile of nb
# Set up a matrix for upper
Fx.CU <- matrix(NA, mat.len, mat.len)
Fy.CU <- matrix(NA, mat.len, mat.len)
Fxy.CU <- matrix(NA, mat.len, mat.len)
Fxx.CU <- matrix(NA, mat.len, mat.len)
Fyy.CU <- matrix(NA, mat.len, mat.len)


# and lower bounds
Fx.CL <- matrix(NA, mat.len, mat.len)
Fy.CL <- matrix(NA, mat.len, mat.len)
Fxy.CL <- matrix(NA, mat.len, mat.len)
Fxx.CL <- matrix(NA, mat.len, mat.len)
Fyy.CL <- matrix(NA, mat.len, mat.len)

for(ri in 1:mat.len) {
  for(cj in 1:mat.len) {
    Fx.CU[ri, cj] <- quantile(grad.Fx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
    Fx.CL[ri, cj] <- quantile(grad.Fx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
    Fy.CU[ri, cj] <- quantile(grad.Fy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
    Fy.CL[ri, cj] <- quantile(grad.Fy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
    Fxy.CU[ri, cj] <- quantile(grad.Fxy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
    Fxy.CL[ri, cj] <- quantile(grad.Fxy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
    Fxx.CU[ri, cj] <- quantile(grad.Fxx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
    Fxx.CL[ri, cj] <- quantile(grad.Fxx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]  
    Fyy.CU[ri, cj] <- quantile(grad.Fyy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
    Fyy.CL[ri, cj] <- quantile(grad.Fyy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]  
  }
  cat(paste("[", ri, ",", cj, "] of ", mat.len, "x", mat.len," \n", sep = ""))
}
#
####################
# Predicted values #
####################
#
## RESPONSE VALUES ##
respgam <- summary(gam1)$formula # best.gam() formula
gam2 <- gam(respgam, data = df.ns,   method = "GCV.Cp", se = T)
#
# xt.ran <- range(df.ns$x)
# yt.ran <- range(df.ns$y)
xt.ran <- range(hull$x)
yt.ran <- range(hull$y)

xtseq <- seq(xt.ran[1], to = xt.ran[2], length.out = mat.len)
ytseq <- seq(yt.ran[1], to = yt.ran[2], length.out = mat.len)

newDt <- expand.grid(list(x = xtseq, y = ytseq))
predt <- predict.gam(gam2, newDt, type = "response")
resp.mat <- matrix(predt, mat.len, mat.len)

## SCALED AND CENTERED VALUES ##
pred.grad <- my.grad(the.gam = gam1, xs = xseq, ys = yseq, eps.mat = eps.list)
newD <- expand.grid(list(x = xseq, y = yseq))
pred <- predict.gam(gam1, newD, type = "response")
pred.mat <- matrix(pred, mat.len, mat.len)

Fx.mat <- matrix(pred.grad$x$grad, mat.len, mat.len)
Fy.mat <- matrix(pred.grad$y$grad, mat.len, mat.len)
Fxy.mat <- matrix(pred.grad$xy$grad, mat.len, mat.len)
Fxx.mat <- matrix(pred.grad$xx$grad, mat.len, mat.len)
Fyy.mat <- matrix(pred.grad$yy$grad, mat.len, mat.len)
#
est.det <- Fxx.mat * Fyy.mat - (Fxy.mat)^2
det.CU <- Fxx.CU * Fyy.CU - (Fxy.CU)^2
det.CL <- Fxx.CL * Fyy.CL - (Fxy.CL)^2
#
save(gam1,                          # Original model
     df, df.ns,                     # S&C and Original data      
     xtseq, ytseq,                  # Original seq
     xseq, yseq,                    # S&C seq   
     grad.Fx, grad.Fy, grad.Fxx, grad.Fyy, grad.Fxy, #
     Fx.CU, Fx.CL,                  #  
     Fy.CU, Fy.CL,                  #
     Fxx.CU,Fxx.CL,                 #    
     Fyy.CU,Fyy.CL,                 # 
     Fxy.CU,Fxy.CL,                 #
     df.parts, df.name,             # Name controls
     pred.mat,                      # S&C Predicted matrix 
     resp.mat,                      # Response predicted matrix 
     det.CU, det.CL,                # Determinant CI
     hull,                          # Convex hull
     file = paste0(df.name, "_DATA_v01.RDATA"))
}

###########
## PLOTS ##
###########
col.mat <- matrix(NA, mat.len, mat.len)
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CL > 0] <- 1 # Local minimum
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CU < 0] <- 2 # Local maximum
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CU < 0] <- 3 # Saddle point

grad.plot <- paste(figure.dir, df.name, "v_001.png", sep="")
png(file = grad.plot, width = 76, height = 80, units = "mm", res = 600)
# layout(matrix(c(1,2), 1,2,  byrow = T),
#        widths = lcm(c(8.45, 8.45)), heights = lcm(c(4.75, 4.75)),
#        respect = F)
nlev <-  10
levs <- pretty(range(pred.mat), nlev)

par(mar=c(2.15,2.15,1.5,0.25), 
#     oma = c(0,0,0,0),
    mgp=c(1,0.5,0), 
    cex= .75, 
    tck= -0.015, 
    family = "sans")
plot(df$x, df$y, type = "n",
     ylab = "",
     xlab = "")
image(xseq,
      yseq,
      pred.mat,
      col = gray.colors(50, start = 0.5, end = 1, gamma = 1, alpha = NULL),
      add = T)
image(xseq,
      yseq,
      col.mat,
      #       col = "white",
      col= c("gray20", "gray30", "gray40"), 
      add = T)
# contour(x = xseq,
#         y = yseq,
#         z = pred.mat,
#         add = T, drawlabels = F)
#
contour(x = xseq,
        y = yseq,
        z = pred.mat,
        zlim = range(pred.mat, finite = T),
        lwd = seq(0, 2, length.out= length(levs)),
        #         col = gray(10:0/ 11),
        col = "black",
        drawlabels = T,
        method = "edge",
        #         labcex = 1,
        #         vfont = c("serif", "bold"),
        vfont = NULL,
        add = T)
legend("topleft", 
       legend = df.parts[1],
       bty = "n")
dev.off()
#
########################################################################################
# If f(x,y)  is a two-dimensional function that has a relative extremum at a point  
# (x0, y0) and has continuous partial derivatives at this point, 
# then fx(x0,y0) = 0 and fy(x0,y0) = 0. The second partial derivatives test classifies 
# the point as a local maximum or relative minimum.

## For the determinant M(x,y) = det(H(x,y)) = 
# if (a, b) is a critical point of f (that is, fx(a, b) = fy(a, b) = 0)
# i.e., fx.CL < 0 & fy.CL < 0 & fx.CU > 0 & fy.CU >0 # We have 95% certainty that we contain zero...

# 1. if M(a,b) > 0 and Fxx(ab) > 0 then (a,b) is a local minimum of f
# 2. if M(a,b) > 0 and Fxx(ab) < 0 then (a,b) is a local maximum of f
# 3. if M(a,b) < 0 and Fxx(ab) then (a,b) is a saddle point of f
# 4. if M(a,b) = 0 then the second derivative test is inconclusive, 
#    and the point (a,b) could be any of a minimum, maximum or saddle point
# Fxx | Fxy
# Fyx | Fyy
#
persp3d(xseq, yseq, Fx.mat, col = col.mat)
persp3d(xseq, yseq, Fx.CL, col = "RED", add = F)
persp3d(xseq, yseq, Fx.CU, col = "GREEN", add = T)
planes3d(0,0,1)
#
persp3d(xseq, yseq, Fy.mat, col = col.mat)
persp3d(xseq, yseq, Fy.CL, col = "RED", add = T)
persp3d(xseq, yseq, Fy.CU, col = "GREEN", add = T)
planes3d(0,0,1)
#
persp3d(xseq, yseq, est.det, col = col.mat)
persp3d(xseq, yseq, det.CL, col = "RED", add = F)
persp3d(xseq, yseq, det.CU, col = "GREEN", add = T)
planes3d(0,0,1)
