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
  setwd("C:/Users/checcomi/Documents/github/newCal/")
}else{
  setwd("/scratch/project_2000994/calibrations/newCal/")
}

devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/newCal/master/Rsrc/functions.r")
devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/newCal/master/Rsrc/settings.r")

##load data
#lastCal <- "chains/calOut_1.rdata"
#newCal <-  "chains/calOut_2.4.rdata"
# load("input/initPrebas.rdata")
load("inputs/init_set1.rdata")
load("inputs/init_set2.rdata")
load("inputs/init_set3.rdata")
load("inputs/init_set4.rdata")
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


### First calibration
# load(lastCal)
# startValue <- calOut$X
# Z <- calOut$Z
startValue <- matrix(NA,3,npar)
startValue[1,] <- runif(npar,parmin,parmax)
startValue[2,] <- runif(npar,parmin,parmax)
startValue[3,] <- runif(npar,parmin,parmax)
startValue[1,1:(length(parSel)*3)] <- c(pCROB[parSel,1],pCROB[parSel,2],pCROB[parSel,3])

settings = list(iterations = iters, startValue = startValue,
                #Z=Z,#burnin=0,message=FALSE)
                message=FALSE)

# ### Create Bayesian Setup
if(vLocal){
  bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                       lower = parmin, upper = parmax,names = parnam
                                       ,parallel = 4)
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
