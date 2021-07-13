# devtools::install_github("checcomi/Rprebasso", force=TRUE,ref="newA")
library(data.table)
library(abind)
library(BayesianTools)
library(plyr)

####set project library if running on CSC
vLocal <- TRUE #### flag for runs on CSC(FALSE) or on laptop(TRUE)
if(!vLocal){
  setwd("/scratch/project_2000994/calibrations/all")
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)
source("")

set1 <- 1:300
set2 <- 301:600
set3 <- 601:922
set4 <- 923:948

# load("/scratch/project_2000994/calibrations/all/inputs/initPrebas.rdata")
load("inputs/initPrebasAll.rdata")

load("outCal/pMAP.rdata")

###initialize parameters
parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)

initPrebas$pCROBAS <- pCROB
initPrebas$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
initPrebas$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
initPrebas$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

####process first set
init_set1 <- subInit(initPrebas, set1)
Hdata_s1 <- subSetData(outdata_H,set1,obs_H)
Ddata_s1 <- subSetData(outdata_D,set1,obs_D)
Bdata_s1 <- subSetData(outdata_B,set1,obs_B)
Vdata_s1 <- subSetData(outdata_V,set1,obs_V)
Hcdata_s1 <- subSetData(outdata_Hc,set1,obs_Hc)

save(init_set1,Hdata_s1,Ddata_s1,Bdata_s1,Vdata_s1,Hcdata_s1,
     file="inputs/init_set1.rdata")


init_set2 <- subInit(initPrebas, set2)
Hdata_s2 <- subSetData(outdata_H,set2,obs_H)
Ddata_s2 <- subSetData(outdata_D,set2,obs_D)
Bdata_s2 <- subSetData(outdata_B,set2,obs_B)
Vdata_s2 <- subSetData(outdata_V,set2,obs_V)
Hcdata_s2 <- subSetData(outdata_Hc,set2,obs_Hc)

save(init_set2,Hdata_s2,Ddata_s2,Bdata_s2,Vdata_s2,Hcdata_s2,
     file="inputs/init_set2.rdata")


init_set3 <- subInit(initPrebas, set3)
Hdata_s3 <- subSetData(outdata_H,set3,obs_H)
Ddata_s3 <- subSetData(outdata_D,set3,obs_D)
Bdata_s3 <- subSetData(outdata_B,set3,obs_B)
Vdata_s3 <- subSetData(outdata_V,set3,obs_V)
Hcdata_s3 <- subSetData(outdata_Hc,set3,obs_Hc)

save(init_set3,Hdata_s3,Ddata_s3,Bdata_s3,Vdata_s3,Hcdata_s3,
     file="inputs/init_set3.rdata")

init_set4 <- subInit(initPrebas, set4)
Hdata_s4 <- subSetData(outdata_H,set4,obs_H)
Ddata_s4 <- subSetData(outdata_D,set4,obs_D)
Bdata_s4 <- subSetData(outdata_B,set4,obs_B)
Vdata_s4 <- subSetData(outdata_V,set4,obs_V)
Hcdata_s4 <- subSetData(outdata_Hc,set4,obs_Hc)

save(init_set4,Hdata_s4,Ddata_s4,Bdata_s4,Vdata_s4,Hcdata_s4,
     file="inputs/init_set4.rdata")
