# devtools::install_github("checcomi/Rprebasso", force=TRUE,ref="newA")
library(data.table)
library(abind)
library(BayesianTools)
library(plyr)

####set project library if running on CSC
vLocal <- FALSE #### flag for runs on CSC(FALSE) or on laptop(TRUE)
if(!vLocal){
  setwd("/scratch/project_2000994/calibrations/all")
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

# function to subset Data
subSetData <- function(outdata,setX,obs){
  selX <- which(outdata[,1] %in% setX)
  if(length(selX)>0){
    obsSub <- obs[selX]
    outdataSub <- outdata[selX,]
    outdataSub[,1] <- outdataSub[,1] - (min(setX)-1)
    return(list(obs=obsSub,outData= outdataSub))
  }
}

####Function to Subset initPrebas 
subInit <- function(initPrebas,setX){
  nYears <- initPrebas$nYears[setX]
  siteInfo=initPrebas$siteInfo[setX,]
  defaultThin=initPrebas$defaultThin[setX]
  ClCut = initPrebas$ClCut[setX]
  multiInitVar = initPrebas$multiInitVar[setX,,]
  yassoRun <- initPrebas$yassoRun[setX]
  multiThin = initPrebas$thinning[setX,,]
  multiNthin = initPrebas$nThinning[setX]
  
  climIDx <- sort(unique(siteInfo[,2]))
  siteInfo[,2] <- match(siteInfo[,2],climIDx)
  maxYears <- max(initPrebas$nYears[setX])
  CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)
  
  for(i in 1:length(climIDx)){
    ij <- climIDx[i]
    PAR[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,1]))
    TAir[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,2]))
    Precip[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,4]))
    VPD[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,3]))
    CO2[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,5]))
  }
  
  init_setX <- InitMultiSite(nYearsMS = nYears,
                             siteInfo=siteInfo,
                             # pCROBAS = pCROBAS,
                             # litterSize = litterSize,#pAWEN = parsAWEN,
                             defaultThin=defaultThin,
                             ClCut = ClCut,
                             multiInitVar = multiInitVar,
                             PAR = PAR,
                             TAir= TAir,
                             VPD= VPD,
                             Precip= Precip,
                             CO2= CO2,
                             yassoRun = yassoRun,#lukeRuns = initPrebas$lukeRuns,
                             # initCLcutRatio = initCLcutRatio
                             multiThin = multiThin,
                             multiNthin = multiNthin
  )
  return(init_setX)
}


set1 <- 1:300
set2 <- 301:600
set3 <- 601:922


# load("/scratch/project_2000994/calibrations/all/inputs/initPrebas.rdata")
load("inputs/initPrebas.rdata")

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

