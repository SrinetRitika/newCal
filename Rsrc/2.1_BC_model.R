# vLocal <- FALSE
if(!vLocal){
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}

library(parallel)
library(BayesianTools)
library(coda)
# library(lhs)

###load data for initialization
if(newV){
  devtools::install_github("ForModLabUHel/Rprebasso", ref="newVersion")
}else{
  devtools::install_github("ForModLabUHel/Rprebasso", ref="master")
}
library(Rprebasso)

if(vLocal){
  setwd("C:/Users/checcomi/Documents/github/newCal/")
}else{
  setwd("/scratch/project_2000994/calibrations/newCal/")
}

devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/newCal/master/Rsrc/functions.r")
devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/newCal/master/Rsrc/settings.r")

# ###load data for initialization
# load(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/init_set1.rdata'))
# load(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/init_set2.rdata'))
# load(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/init_set3.rdata'))
# load(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/init_set4.rdata'))
if(newV){
  load('inputs/init_set1_newVersion.rdata')
  load('inputs/init_set2_newVersion.rdata')
  load('inputs/init_set3_newVersion.rdata')
  load('inputs/init_set4_newVersion.rdata')
  load('inputs/init_set5Flux_newVersion.rdata')
}else{
  load('inputs/init_set1.rdata')
  load('inputs/init_set2.rdata')
  load('inputs/init_set3.rdata')
  load('inputs/init_set4.rdata')
  load('inputs/init_set5Flux.rdata')
}

vapu_S<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_spruce.csv'))
nData_S <- length(vapu_S$plotNo)
vapu_P<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_pine.csv'))
nData_P <- length(vapu_P$plot)



# ### Create Bayesian Setup
if(vLocal){
  bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                       lower = parmin, upper = parmax,names = parnam,
                                       parallel = nworkers)
}else{
  bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                       lower = parmin, upper = parmax,
                                       names = parnam)
}

### First calibration
if(calSet == 0){
  startValue <- matrix(NA,3,npar)
  startValue[1,] <- runif(npar,parmin,parmax)
  startValue[2,] <- runif(npar,parmin,parmax)
  startValue[3,] <- runif(npar,parmin,parmax)
  startValue[1,1:(length(parSel)*3)] <- c(pCROB[parSel,1],pCROB[parSel,2],pCROB[parSel,3])
  
  settings = list(iterations = iters, f=fX, 
                  startValue = startValue,message=FALSE)
  
  calOut <- runMCMC(bayesianSetup = bayesianSetup, 
                    sampler = calAlg, settings = settings)
}else{
  ####continue calibration
  load(lastCal)
  startValue <- calOut$X
  
  settings = list(iterations = iters, startValue = startValue,
                  burnin=0,message=FALSE)
  
  calOut <- runMCMC(bayesianSetup = bayesianSetup,
                    sampler = calAlg, settings = settings)
}
save(calOut,file=newCal)