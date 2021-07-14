

####Run model in parallel
library(ggplot2)
# library(Metrics)
library(data.table)
library(parallel)
library(BayesianTools)
library(ggpubr)

####set project library if running on CSC
vLocal <- TRUE #### flag for runs on CSC(FALSE) or on laptop(TRUE)
if(!vLocal){
  setwd("/scratch/project_2000994/calibrations/all")
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

# setwd("C:/Users/minunno/Documents/research/calibrationNew/all")

source("settings.r") ###run using url
source("functions.r") ###run using url

# load("inputs/initPrebas_old growth.rdata")
load("inputs/init_set1.rdata")
load("inputs/init_set2.rdata")
load("inputs/init_set3.rdata")
load("inputs/init_set4.rdata")

load("outCal/pMAP.rdata")



###change parameters
###init_set1
init_set1$pCROBAS <- pCROB
init_set1$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set1$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set1$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
###init_set2
init_set2$pCROBAS <- pCROB
init_set2$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set2$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set2$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
###init_set3
init_set3$pCROBAS <- pCROB
init_set3$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set3$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set3$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
###init_set4
init_set4$pCROBAS <- pCROB
init_set4$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set4$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set4$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

###runModel
modOut1 <- multiPrebas(init_set1)
modOut2 <- multiPrebas(init_set2)
modOut3 <- multiPrebas(init_set3)
modOut4 <- multiPrebas(init_set4)

###calculate likelihoods
likelihood1(pMAP)
likelihood2(pMAP)
likelihood3(pMAP)
likelihood4(pMAP)


 nCores <- 4
 sets <- 1:4
 likelihoods <- list()
 likelihoods[[1]] <- likelihood1
 likelihoods[[2]] <- likelihood2
 likelihoods[[3]] <- likelihood3
 likelihoods[[4]] <- likelihood4
 modOut <- list()
# 
###Run model using 1 initialization object
startX <- Sys.time()
modOutAll <- multiPrebas(initPrebas)
endX <- Sys.time()
timeAll = endX- startX


library("parallel")
startX <- Sys.time()
logLike <- mclapply(sets, function(jx) {
  likelihoods[[jx]](pMAP)  ## Do nothing for 10 seconds
})#, mc.cores = nCores)  
endX <- Sys.time()
timeX = endX- startX


startX <- Sys.time()
likelihood1(pMAP)
likelihood2(pMAP)
likelihood3(pMAP)
likelihood4(pMAP)
endX <- Sys.time()
timeX = endX- startX




