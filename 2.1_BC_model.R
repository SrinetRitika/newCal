vLocal <- FALSE
if(!vLocal){
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}

library(parallel)
library(BayesianTools)
library(coda)
# library(lhs)

devtools::install_github("ForModLabUHel/Rprebasso", ref="master")
library(Rprebasso)

if(vLocal){
  setwd("C:/Users/minunno/Documents/research/calibrationNew/all")
}else{
  setwd("/scratch/project_2000994/calibrations/all")
}
  
##load data
#lastCal <- "chains/calOut_1.rdata"
#newCal <-  "chains/calOut_2.4.rdata"
load("input/initPrebas.rdata")

# if(vLocal){
#   nworkers <- 4 # required cores
#   iters <- 1.2e5
# }else{
#   nworkers <- 30 # required cores
#   iters <- 6e2
#   #Make a cluster with nworkers
#   cl <- parallel::makeCluster(nworkers)
# }



##settings
# nCores <- 30
# fileNamelastCal <- "chains/resDE_n0.RData"
# fileNameSave <- "chains/resDE_n0.1.RData"

# load("input/pValues.rdata")
###process parameters
param_all <- read.csv('input/parameters.csv')
parSel <- c(1:18,31:34,38,41)
parmod <- param_all[,2] 
parmin <- param_all[,3] 
parmax <- param_all[,4] 
parnam <- as.character(param_all[,1])
npar <- length(parmod)

##load Start pValues matrix
# load("input/startValue.rdata")
nparCROB <- length(parSel)
outdata_B2 <- outdata_B;outdata_B2[,5] <- 2
outdata_V2 <- outdata_V;outdata_V2[,5] <- 2
initPrebas$pCROBAS <- pCROB
# initPrebas$pCROBAS[12,] <- initPrebas$pCROBAS[12,]-1

# Heavy tailed noraml distribution
Sivia_log<-function(diff,sd){
  # sd[which(sd<=0)]<-11e-6e-6
  diff[which(abs(diff)<=1e-6)]<-1e-6
  R2<-(diff/sd)^2
  prob<-1/(sd*(pi*2)^0.5)*(1-exp(-R2/2))/R2
  log(prob)
}

likelihood <- function(pValues){
  nparCROB=24
  initPrebas$pCROBAS[parSel,1] <- pValues[1:nparCROB]
  initPrebas$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
  initPrebas$pCROBAS[parSel,3] <- pValues[(nparCROB*2+1):(nparCROB*3)]
  
  output <- multiPrebas(initPrebas)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {

  out_B <-  output[outdata_B] + output[outdata_B2]
  out_V <-  output[outdata_V] + output[outdata_V2]
  diff_H <- output[outdata_H]-obs_H
  diff_Hc <- output[outdata_Hc]-obs_Hc
  diff_D <- output[outdata_D]-obs_D
  diff_B <- out_B-obs_B
  diff_V <- out_V-obs_V
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB+1)]+pValues[(nparCROB+2)]*output[outdata_H]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB+3)]+pValues[(nparCROB+4)]*output[outdata_D]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB+5)]+pValues[(nparCROB+6)]*output[outdata_B]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB+7)]+pValues[(nparCROB+8)]*output[outdata_Hc]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB+9)]+pValues[(nparCROB+10)]*output[outdata_V]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))

  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }

  return(loglikelihood)
}


### First calibration
load(lastCal)
startValue <- calOut$X
Z <- calOut$Z
# startValue <- matrix(NA,3,npar)
# startValue[1,] <- runif(npar,parmin,parmax)
# startValue[2,] <- runif(npar,parmin,parmax)
# startValue[3,] <- runif(npar,parmin,parmax)
# startValue[1,1:72] <- c(pCROB[parSel,1],pCROB[parSel,2],pCROB[parSel,3])

settings = list(iterations = iters, startValue = startValue,
                Z=Z,#burnin=0,message=FALSE)
                message=FALSE)

# ### Create Bayesian Setup
if(vLocal){
  bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                       lower = parmin, upper = parmax,names = parnam
                                       ,parallel = nworkers)
}else{
  bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                       lower = parmin, upper = parmax,
                                       names = parnam)#,
                                       # parallel = "external")
}
 

calOut <- runMCMC(bayesianSetup = bayesianSetup, 
                  sampler = "DREAMzs", settings = settings)
save(calOut,file=newCal)

# ###Some plot
# pMAP <- MAP(calOut)[[1]]

# source("runModel.r")

# pdf("results.pdf")
# plot(calOut)
# correlationPlot(calOut)
# marginalPlot(calOut, prior = TRUE)
# gelmanDiagnostics(calOut, plot = T)
# dev.off()



# 
# 
#  # settings = list(iterations = iters, startValue = startValue)
#  #  oldw <- getOption("warn")
#  # options(warn = -1)
#  #stX <- Sys.time()
#  #resDE <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
#  # options(warn = oldw)
#  # stX <- Sys.time()
#  # 
#  # save(resDE, file = fileNameSave)
#  # enX <- Sys.time()
#  # timeRun <- enX - stX
#  # print(timeRun)
# 
#  
#  
# # ####continue calibration for parallel supercomp
#  load(fileNamelastCal)
# 
#  
# # X <- matrix(NA,3,82)
# # X[,1:23] <- resDE$X[,1:23]
# # X[,c(24,48,72)] <- 1
# # X[,25:47] <- resDE$X[,24:46]
# # X[,49:71] <- resDE$X[,47:69]
# # X[,73:82] <- resDE$X[,70:79]
#  startValue <- resDE$X
#  
#  Z <- resDE$Z
#  settings = list(iterations = iters, startValue = startValue,
#                  # Z=Z,burnin=0,message=FALSE)
#                   burnin=0,message=TRUE)
# 
#  stX <- Sys.time()
# 
#  resDE <- runMCMC(bayesianSetup = bayesianSetup, 
#                   sampler = "DREAMzs", settings = settings)
# 
#  save(resDE, file = fileNameSave) #save(resDE, file = "resDEtest.rdata")
#  enX <- Sys.time()
#  timeRun <- enX - stX
#  print(timeRun)
# 
# #Stop the cluster
# if(!vLocal){
#    stopCluster(cl)
#  }
