

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


# load("inputs/initPrebas_old growth.rdata")
load("inputs/initPrebas.rdata")
load("inputs/init_set1.rdata")
load("inputs/init_set2.rdata")
load("inputs/init_set3.rdata")

load("outCal/pMAP.rdata")


parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)

initPrebas$pCROBAS <- pCROB
initPrebas$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
initPrebas$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
initPrebas$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

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


nCores <- 3
runs <- 1:3
initX <- list()
initX[[1]] <- init_set1
initX[[2]] <- init_set2
initX[[3]] <- init_set3
modOut <- list()

###Run model using 1 initialization object
startX <- Sys.time()
modOutAll <- multiPrebas(initPrebas)  
endX <- Sys.time()
timeAll = endX- startX

####Run model for splitted sites
modOut1 <- multiPrebas(init_set1)  
modOut2 <- multiPrebas(init_set2)  
modOut3 <- multiPrebas(init_set3)  

####extract and compare data
extraxtData <- function(modOutAll,modOut1,modOut2,modOut3,
                        outdata,data_s1,data_s2,data_s3,obs_all){
  sim_all <- modOutAll$multiOut[outdata]
  sim_set1 <- modOut1$multiOut[data_s1$outData]
  sim_set2 <- modOut2$multiOut[data_s2$outData]
  sim_set3 <- modOut3$multiOut[data_s3$outData]
  
  plot(sim_all,obs_all)
  points(sim_set1,data_s1$obs,pch=20,col=2)
  points(sim_set2,data_s2$obs,pch=20,col=3)
  points(sim_set3,data_s3$obs,pch=20,col=4)
  
  plot(obs_all,c(data_s1$obs,data_s2$obs,data_s3$obs))
  plot(sim_all,c(sim_set1,sim_set2,sim_set3))
}

extraxtData(modOutAll,modOut1,modOut2,modOut3,
            outdata_B,Bdata_s1,Bdata_s2,Bdata_s3,obs_B)
# 
# 
# ###Run in parallel only working on linux
# print("I'm here")
# startX <- Sys.time()
# modOutParl <- mclapply(runs, function(jx) {
#   modOut[[jx]] <- multiPrebas(initX[[jx]])  
# }, mc.cores = nCores)      
# endX <- Sys.time()
# timeParl = endX- startX
# # print(timeX)
# 
# print("I'm here 2")
# startX <- Sys.time()
# modOutAll <- multiPrebas(initPrebas)  
# endX <- Sys.time()
# timeAll = endX- startX
# # print(timeX)
# 
# print("I'm here 3")
# startX <- Sys.time()
# modOut1 <- multiPrebas(init_set1)  
# modOut2 <- multiPrebas(init_set2)  
# modOut3 <- multiPrebas(init_set3)  
# endX <- Sys.time()
# timeSerial = endX- startX
# 
# print(paste("parallel", timeParl))
# print(paste("serial", timeSerial))
# print(paste("All", timeAll))
# 
