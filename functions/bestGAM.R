################
# FORMULA INFO #
################
# PURPOSE: Function that uses BIC model selection to identify the "best" model for a data series
# and creates a data.frame with with all relevant model outputs.
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# VERSION: 0.1
#

# formulaList is a list(c("F1"...)) of GAM models (where, F1 <- formula(z ~ s(x,  k = k.2))...)
# dataFrame is the data.frame to be analyzed, where df <- data.frame("x" = x, y = y, "z" = z)
# gamMethod (see ?mgcv for description of methods) "GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", and "P-ML"

best.gam <- function(formulaList, dataFrame, gamMethod) {
# Make sure "mgcv" is loaded:
  if(require("mgcv", character.only = T)) {
  } else {
    print("trying to install mgcv")
    install.packages("mgcv")
    if(require(mgcv)) {
      print("mgcv installed and loaded")
    } else {
      stop("could not install mgcv")
    }
  }
  gam.fits <- data.frame()
  for(f in formulaList){
    gam.f <- gam(get(f), method = gamMethod, se = T, data = dataFrame)
## Collect appropriate P-values and estimated degrees of freedom for each term
 # For models with 2 smoothing terms
    if(length(summary(gam.f)$s.pv) == 2) { 
      X_pvalue <- round(summary(gam.f)$s.pv, 3)[1]
      Y_pvalue <- round(summary(gam.f)$s.pv, 3)[2]
      X_EDF    <- round(summary(gam.f)$edf, 3)[1]
      Y_EDF    <- round(summary(gam.f)$edf, 3)[2]
    } 
 # For models with linear and nonlinear components
    if(length(summary(gam.f)$s.pv) == 1 & length(summary(gam.f)$pTerms.pv) == 1) { 
      # Linear "x" and nonlinear "y"
      if(all(dimnames(summary(gam.f)$pTerms.pv) == "x")) {
        X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
        Y_pvalue <- round(summary(gam.f)$s.pv, 3)
        X_EDF    <- NA
        Y_EDF    <- round(summary(gam.f)$edf, 3)
      }
      # Nonlinear "x" and linear "y"
      if(all(dimnames(summary(gam.f)$pTerms.pv) == "y")) {
        X_pvalue <- round(summary(gam.f)$s.pv, 3)
        Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
        X_EDF    <- round(summary(gam.f)$edf, 3)
        Y_EDF    <- NA
      } 
    }
 # For models with linear "x" and "y"
    if(length(summary(gam.f)$pTerms.pv) == 2) {
      X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)[1]
      Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)[2]
      X_EDF    <- NA
      Y_EDF    <- NA
    }
 # For models with a single driver: s(x), s(y), x, and y
    if(length(summary(gam.f)$s.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "x")) {
      X_pvalue <- round(summary(gam.f)$s.pv, 3)
      Y_pvalue <- NA
      X_EDF    <- round(summary(gam.f)$edf, 3)
      Y_EDF    <- NA
    }
    if(length(summary(gam.f)$s.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "y")) {
      X_pvalue <- NA
      Y_pvalue <- round(summary(gam.f)$s.pv, 3)
      X_EDF    <- NA
      Y_EDF    <- round(summary(gam.f)$edf, 3)
    }
    if(length(summary(gam.f)$pTerms.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "x")) {
      X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
      Y_pvalue <- NA
      X_EDF    <- NA
      Y_EDF    <- NA
    }
    if(length(summary(gam.f)$pTerms.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "y")) {
      X_pvalue <- NA
      Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
      X_EDF    <- NA
      Y_EDF    <- NA
    }
## Create a data frame of all GAM for formulaList
    gam.fits <- rbind(gam.fits,
                      data.frame("MODEL" = f,
                                 "AICc"= AICc(gam.f),
                                 "BIC" = BIC(gam.f),
                                 "GCF" = summary(gam.f)$sp.criterion,
                                 "DEV" = summary(gam.f)$dev.expl,
                                 "X_p-value" = X_pvalue, 
                                 "Y_p-value" = Y_pvalue,
                                 "X_EDF" = X_EDF,
                                 "Y_EDF" = Y_EDF,
                                 stringsAsFactors=FALSE))
    assign(paste(f, "gam", sep = "."), gam.f)
  }
return(gam.fits)
}

# Calculate AICc, which corrects for finite sample size Burnham & Anderson 2002
AICc <- function(mod) {
  K.c <- mod$rank
  N.c <- length(mod$residuals)
  AIC.c <- round(mod$aic + (2*K.c*(K.c+1)/(N.c-K.c-1)),3)
  return(AIC.c)
}