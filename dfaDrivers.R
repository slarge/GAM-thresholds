rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: MARSS DFA approach on ecological indicators
# AUTHOR: Scott Large 2012 (revision 2015)
# REVIEWED BY:
# VERSION: 0.3
#

######################
# CHANGES/ ADDITIONS #
######################
# Added:
#
# Done:
#
############
# PACKAGES #
############
library(MARSS)
library(reshape2)
#
##################
# Set up folders #
##################
data.dir <- "data/"
figure.dir <- "output/"
source.dir <- "functions/"
#
#############
# Load data #
#############
set.seed(627)
cc <- read.csv("data/CCIEA-RPW.csv")
# Subset area... initally for coastwide
ccALL <- cc[cc$Coastwide == 1,]
ccALL$year <- as.numeric(ccALL$year)
#
ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
# Wide df with columns as variables
dat.full <- dcast(ccALL, year ~ timeseries, value.var = "value")
#
ind.name <- c("GF_spp_richness_coastwide", "GF-Simp_coastwide", "GF_MTL_coastwide", "Scav_ratio")
driver <- dat.full[!colnames(dat.full) %in% ind.name]
year <- driver$year
#
# transpose data, get rid of YEAR and LANDINGS
df <- driver[!colnames(driver) %in% c("year")]
dat <- t(df)

N.ts <- dim(dat)[1]
# get length of time series
TT <- dim(dat)[2] 

# Take the z-scores of data to normalize
Sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
Mean <- apply(dat,1,mean,na.rm=TRUE)
dat.z <- (dat-Mean)*(1/Sigma)
rownames(dat.z) <- colnames(df)
# save(dri.data, file = "dfaDATA_DRI_v0001.Rdata")

##############
# BASE MODEL #
##############
# set new control params
cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
# set up forms of R matrices
levels.R <- c("diagonal and equal",
              "diagonal and unequal")
model.data = data.frame()

# NOTE: this will take a long time to run! For the purpose of the dfaGAM we will
# start with N.ts -1 = 4
N.ts = 5
for(R in levels.R) {
#     for(m in 1:(N.ts-1)) {
    for(m in 1:(N.ts-1)) {    
        dfa.model = list(A="zero", R=R, m=m)
        kemz = MARSS(dat.z, 
                     model=dfa.model, 
                     control=cntl.list,
                     silent = T,
                     form="dfa",
                     z.score=TRUE)
        model.data = rbind(model.data,
                            data.frame(R=R,
                                       m=m,
                                       logLik=kemz$logLik,
                                       K=kemz$num.params,
                                       AICc=kemz$AICc,                                       
                                       stringsAsFactors=FALSE))
        assign(paste("kemz", m, R, sep="."), kemz)
        cat("Just finished",m,"hidden trend(s) with a",R,"covariance matrix with no covariate \n") # end m loop
        } # end m loop
    } # end R loop
save(file="DFA_drivers.Rdata",list = c("model.data", ls(pattern="^kemz."))) # end R loop
# load("DFA_drivers.Rdata")


##############################################################
# Create table of different structures and numbers of trends #
##############################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)
#
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]
#
best.model <- model.tbl[1,]
fitname = paste("kemz",best.model$m,best.model$R,sep=".")
best.fit = get(fitname)
#
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat
# H.inv = promax(coef(best.fit, type="matrix")$Z)$rotmat
#
# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv   
# rotate trends
trends.rot = solve(H.inv) %*% best.fit$states
ts.trends = cbind(dat.full$year, t(trends.rot))
colnames(ts.trends) <- c("YEAR", "T1", "T2", "T3")
#
save(driver, dat, dat.z, ts.trends, file = paste0(data.dir, "dfaTrends_v01.RDATA"))
# write.csv(model.tbl, paste0(figure.dir, "dfaResults_v0002.csv"), row.names = F)