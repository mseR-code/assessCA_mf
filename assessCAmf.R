## ------------------------------------------------------------------
## Control script for multifleet assessCA TMB
## model.
##
##
## Author: SDN Johnson
## Date: Oct 22, 2018
##
## Date last revised:
##
## Reads in data for a stock that has multiple commercial fleets and
## fishery (in)dependent surveys, and runs the assessCAmf TMB model.
## Designed to test a new model for mseR, and potentially for other
## OM conditioning applications.
##
## ToDo:
##   1. Control file capability
##   2. Separate simulation model that generates data for
##      this fitting script - not super necessary but req'd
##      for later expansion to OM conditioning
##
## ------------------------------------------------------------------

# Handling script for the TMB model assessDD.cpp
# Source tools, load packages
source("mseRtools.r")
source("assessCAfuns.R")
library(TMB)            
library(RColorBrewer)        

# Compile and load TMB model
compile("assessCA.cpp")               # Compile the C++ file
dyn.load(dynlib("assessCA"))          # Dynamically link the C++ code


# Load the data
ages  <- lisread("opModAges.dat")
catch <- lisread("opModCatch.dat")
idx   <- lisread("opModIndex.dat")

nT    <- 54

# Put data into correct shape
C_tg   <- catch$landCatchMatrix[1:nT,3:7]/1e3
I_tg   <- matrix(-1, nrow = nT, ncol = 5)
I_tg[,c(1,4,5)] <- t(idx$idxSeries)

# Ages
A_atg <- array(-1, dim = c(35,nT,5))
A_atg[2:35,,1] <- t(ages$ageobsProp_f1)
A_atg[2:35,,4] <- t(ages$ageobsProp_f2)
A_atg[2:35,,5] <- t(ages$ageobsProp_f3)

for(g in 1:5)
  for(t in 1:nT)
    if(abs(sum(A_atg[,t,g])) < 1e-1 ) 
      A_atg[1,t,g] <- 0

firstRecDev <- 10

dat <- list(  I_tg = I_tg,
              C_tg = C_tg,
              A_atg = A_atg,
              survType_g = c(1,-1,-1,1,1),
              indexType_g = c(0,-1,-1,0,0),
              calcIndex_g = c(1,0,0,1,1),
              selType_g = c(0,0,0,0,0),
              fleetTiming = c(.49,.5,.51,.74,.76),
              initCode = c(0),
              posPenFactor = c(1e3),
              firstRecDev = firstRecDev )

par <- list(  lnB0 = 4,
              logit_ySteepness = 0,
              lnM = -2.7,
              log_initN_mult = rep(0,35),
              lnSelAlpha_g = rep(1.6,5),
              lnSelBeta_g = rep(2.4,5),
              lntauAge_g = rep(-1.6,5),
              effSampleSize_g = rep(100,5),
              recDevs_t = rep(0,nT-firstRecDev + 1 ),
              lnsigmaR = 0,
              omegaM_t = rep(0,nT-1),
              lnsigmaM = -1,
              agetau2IGa = rep(1,5),
              agetau2IGb = rep(0.04,5),
              obstau2IGa = rep(1,5),
              obstau2IGb = rep(0.04,5),
              sig2RPrior = c(1,2),
              sig2MPrior = c(1,0.04),
              rSteepBetaPrior = c(.55,.005),
              initMPrior = c(.06,.006),
              aMat = c(5,12),
              Linf = c(70),
              L1 = c(32.5),
              vonK = c(.275),
              lenWt = c(1.04e-5,3.07))


map <- list(  # lnB0 =factor(NA),
              # logit_ySteepness = factor(NA),
              # lnM = factor(NA),
              log_initN_mult = factor(rep(NA,35)),
              #lnSelAlpha_g = factor(rep(NA,5)),
              #lnSelBeta_g = factor(rep(NA,5)),
              #recDevs_t = factor(rep(NA,nT-firstRecDev + 1)),
              lnsigmaR = factor(NA),
              omegaM_t = factor(rep(NA,nT-1)),
              lnsigmaM = factor(NA),
              agetau2IGa = factor(rep(NA,5)),
              agetau2IGb = factor(rep(NA,5)),
              obstau2IGa = factor(rep(NA,5)),
              obstau2IGb = factor(rep(NA,5)),
              sig2RPrior = factor(c(NA,NA)),
              sig2MPrior = factor(c(NA,NA)),
              rSteepBetaPrior = factor(c(NA,NA)),
              initMPrior = factor(c(NA,NA)),
              aMat = factor(c(NA,NA)),
              Linf = factor(c(NA)),
              L1 = factor(c(NA)),
              vonK = factor(c(NA)),
              lenWt = factor(c(NA,NA)) )



ctrl <- list( eval.max = 1000, iter.max = 1000 )


objFE <- MakeADFun( data=dat,parameters=par,map=map, 
                    random = NULL, silent = TRUE )
                    # random=c("omegaR_t") )

fitFE <- try( nlminb (  start     = objFE$par,
                        objective = objFE$fn,
                        gradient  = objFE$gr,
                        control   = ctrl ) )

repFE <- objFE$report()
sdrepFR <- summary(sdreport(objFE))

