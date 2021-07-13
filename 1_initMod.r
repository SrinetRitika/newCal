# devtools::install_github("checcomi/Rprebasso", force=TRUE,ref="newA")
library(data.table)
# library(Rprebasso)
library(abind)
library(BayesianTools)
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

load("C:/Users/manhu/OneDrive/Research/PREBAS simulation/newCal-master/Pine/initPrebas.rdata")
pine <- initPrebas
pine_B <- outdata_B
pine_D <- outdata_D
pine_H <- outdata_H
pine_V <- outdata_V
obs_pB <- obs_B
obs_pD <- obs_D
obs_pH <- obs_H
obs_pV <- obs_V
load("C:/Users/manhu/OneDrive/Research/PREBAS simulation/newCal-master/Spruce/initPrebas.rdata")
spruce <- initPrebas
spruce_B <- outdata_B
spruce_D <- outdata_D
spruce_H <- outdata_H
spruce_Hc <- outdata_Hc
spruce_V <- outdata_V
obs_sB <- obs_B
obs_sD <- obs_D
obs_sH <- obs_H
obs_sHc <- obs_Hc
obs_sV <- obs_V
load("C:/Users/manhu/OneDrive/Research/PREBAS simulation/newCal-master/Birch/initPrebas.rdata")
nfi <- initPrebas
nfi_B <- outdata_B
nfi_D <- outdata_D
nfi_H <- outdata_H
obs_nB <- obs_B
obs_nD <- obs_D
obs_nH <- obs_H
load("C:/Users/manhu/OneDrive/Research/PREBAS simulation/newCal-master/inputs/initPrebas_old growth.rdata")
old <- initPrebas
old_B <- as.matrix(outdata_B)
old_D <- as.matrix(outdata_D)
old_H <- as.matrix(outdata_H)
old_Hc<- as.matrix(outdata_Hc)
old_V <- as.matrix(outdata_V)
obs_oB <- obs_B
obs_oD <- obs_D
obs_oH <- obs_H
obs_oHc <- obs_Hc
obs_oV <- obs_V

initPrebas$nSites <- pine$nSites + spruce$nSites + nfi$nSites
initPrebas$multiOut <- array(0,dim=c(initPrebas$nSites,
                                     max(pine$nYears,spruce$nYears,nfi$nYears),
                                     46,3,2))
initPrebas$nClimID <- pine$nClimID + spruce$nClimID + nfi$nClimID
# initPrebas$nLayers
initPrebas$maxYears <- max(pine$nYears,spruce$nYears,nfi$nYears)
initPrebas$maxThin <- max(pine$maxThin,spruce$maxThin,nfi$maxThin)
initPrebas$nYears <- c(pine$nYears,spruce$nYears,nfi$nYears)
initPrebas$areas <- rep(1,initPrebas$nSites)

spruce$thinning[,,2] <- spruce$thinning[,,2] + pine$nSites
initPrebas$thinning <- array(0, dim=c(initPrebas$nSites,initPrebas$maxThin,9))
initPrebas$thinning[,,9] <- -999
initPrebas$thinning[1:pine$nSites,,1:8] <- pine$thinning
initPrebas$thinning[(pine$nSites+1):(pine$nSites+spruce$nSites),1:spruce$maxThin,1:8] <- spruce$thinning

nfi$siteInfo[,2] <- nfi$siteInfo[,2] + pine$nClimID + spruce$nClimID
spruce$siteInfo[,2] <- spruce$siteInfo[,2] + pine$nClimID
nfi$siteInfo[,1] <- nfi$siteInfo[,1] + pine$nSites + spruce$nSites
spruce$siteInfo[,1] <- spruce$siteInfo[,1] + pine$nSites

initPrebas$siteInfo <- rbind(pine$siteInfo,spruce$siteInfo,nfi$siteInfo)

initPrebas$nThinning <- c(pine$nThinning,spruce$nThinning,nfi$nThinning)

initPrebas$fAPAR <- matrix(0.7,initPrebas$nSites,initPrebas$maxYears)
initPrebas$initClearcut <- rbind(pine$initClearcut,spruce$initClearcut,nfi$initClearcut)
initPrebas$fixBAinitClarcut <- rep(1,initPrebas$nSites)
initPrebas$initCLcutRatio <- matrix(0,initPrebas$nSites,3)

initPrebas$ETSy <- matrix(0,initPrebas$nClimID,initPrebas$maxYears)
initPrebas$ETSy[1:pine$nClimID,] <- pine$ETSy
initPrebas$ETSy[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears] <- spruce$ETSy
initPrebas$ETSy[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears] <- nfi$ETSy

initPrebas$P0y <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,2))
initPrebas$P0y[1:pine$nClimID,,] <- pine$P0y
initPrebas$P0y[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,] <- spruce$P0y
initPrebas$P0y[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,] <- nfi$P0y

initPrebas$multiInitVar <- array(0,dim=c(initPrebas$nSites,7,3))
initPrebas$multiInitVar[1:pine$nSites,,1] <- pine$multiInitVar
initPrebas$multiInitVar[(pine$nSites+1):(pine$nSites + spruce$nSites),,1] <- spruce$multiInitVar
initPrebas$multiInitVar[(pine$nSites+1+ spruce$nSites):(nfi$nSites + pine$nSites + spruce$nSites),,] <- nfi$multiInitVar

initPrebas$weather <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,365,5))
initPrebas$weather[1:pine$nClimID,,,] <- pine$weather
initPrebas$weather[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,,] <- spruce$weather
initPrebas$weather[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,,] <- nfi$weather

initPrebas$soilC <- array(0,dim=c(initPrebas$nSites,initPrebas$maxYears,5,3,3))

initPrebas$weatherYasso <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,3))
# initPrebas$weatherYasso[1:pine$nClimID,,] <- pine$weatherYasso
# initPrebas$weatherYasso[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,] <- spruce$weatherYasso
# initPrebas$weatherYasso[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,] <- nfi$weatherYasso

initPrebas$soilCtot <- matrix(0,initPrebas$nSites,initPrebas$maxYears)

initPrebas$defaultThin <- rep(0,initPrebas$nSites)

initPrebas$ClCut <- rep(0,initPrebas$nSites)

initPrebas$inDclct <- matrix(0,initPrebas$nSites,3)
initPrebas$inAclct <- matrix(0,initPrebas$nSites,3)

initPrebas$dailyPRELES <- array(0,dim=c(initPrebas$nSites,(initPrebas$maxYears*365),3))

initPrebas$yassoRun <- rep(0,initPrebas$nSites)

initPrebas$PREBASversion <- rep(0,initPrebas$nSites)


spruce_B[,1] <- spruce_B[,1] + pine$nSites
spruce_D[,1] <- spruce_D[,1] + pine$nSites
spruce_H[,1] <- spruce_H[,1] + pine$nSites
spruce_Hc[,1] <- spruce_Hc[,1] + pine$nSites
spruce_V[,1] <- spruce_V[,1] + pine$nSites
nfi_B[,1] <- nfi_B[,1] + pine$nSites + spruce$nSites
nfi_D[,1] <- nfi_D[,1] + pine$nSites + spruce$nSites
nfi_H[,1] <- nfi_H[,1] + pine$nSites + spruce$nSites

outdata_B <- rbind(pine_B,spruce_B,nfi_B)
outdata_D <- rbind(pine_D,spruce_D,nfi_D)
outdata_H <- rbind(pine_H,spruce_H,nfi_H)
outdata_Hc <- spruce_Hc
outdata_V <- rbind(pine_V,spruce_V)


obs_B <- c(obs_pB,obs_sB,obs_nB)
obs_D <- c(obs_pD,obs_sD,obs_nD)
obs_H <- c(obs_pH,obs_sH,obs_nH)
obs_Hc <- obs_Hc
obs_V <- c(obs_pV,obs_sV)


CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,initPrebas$nClimID,initPrebas$maxYears*365)

for(i in 1:initPrebas$nClimID){
  PAR[i,] <- as.vector(t(initPrebas$weather[i,,,1]))
  TAir[i,] <- as.vector(t(initPrebas$weather[i,,,2]))
  Precip[i,] <- as.vector(t(initPrebas$weather[i,,,4]))
  VPD[i,] <- as.vector(t(initPrebas$weather[i,,,3]))
  CO2[i,] <- as.vector(t(initPrebas$weather[i,,,5]))
}

multiThin <- initPrebas$thinning
nThinning <- initPrebas$nThinning
initPrebas$siteInfo <- cbind(initPrebas$siteInfo,matrix(pPREL[1:3],nrow(initPrebas$siteInfo),3,byrow = T))
pCROBAS <- pCROB
pCROBAS[12,] <- pCROBAS[12,] - 1

initPrebas <- InitMultiSite(nYearsMS = initPrebas$nYears,
                            siteInfo=as.matrix(initPrebas$siteInfo),
                            pCROBAS = pCROBAS,
                            # litterSize = litterSize,#pAWEN = parsAWEN,
                            defaultThin=initPrebas$defaultThin,
                            ClCut = as.integer(initPrebas$ClCut),
                            multiInitVar = as.array(initPrebas$multiInitVar),
                            PAR = as.matrix(PAR),
                            TAir= as.matrix(TAir),
                            VPD= as.matrix(VPD),
                            Precip= as.matrix(Precip),
                            CO2= as.matrix(CO2),
                            yassoRun = initPrebas$yassoRun,#lukeRuns = initPrebas$lukeRuns,
                            # initCLcutRatio = initCLcutRatio
                            multiThin = as.array(multiThin),
                            multiNthin = as.matrix(nThinning)
)

save(initPrebas,obs_H,obs_D,obs_B,obs_V,obs_Hc,
     outdata_H,outdata_D,outdata_B,outdata_V,outdata_Hc,
     file="input/initPrebas.rdata")


#########together with old growth
nvar<-length(varNames)
initPrebas$nSites <- pine$nSites + spruce$nSites + nfi$nSites+old$nSites
initPrebas$multiOut <- array(0,dim=c(initPrebas$nSites,
                                     max(pine$nYears,spruce$nYears,nfi$nYears,old$nYears),
                                     nvar,3,2))  ## what is 46, 3, 2, 46 should be the number of variables in the multiOut, 3max mum of layers, 2 related to the thinning
initPrebas$nClimID <- pine$nClimID + spruce$nClimID + nfi$nClimID + old$nClimID

# initPrebas$nLayers
initPrebas$maxYears <- max(pine$nYears,spruce$nYears,nfi$nYears,old$nYears)
initPrebas$maxThin <- max(pine$maxThin,spruce$maxThin,nfi$maxThin,old$maxThin)
initPrebas$nYears <- c(pine$nYears,spruce$nYears,nfi$nYears,old$nYears)
initPrebas$areas <- rep(1,initPrebas$nSites)   #what does the area stand for here?

spruce$thinning[,,2] <- spruce$thinning[,,2] + pine$nSites
initPrebas$thinning <- array(0, dim=c(initPrebas$nSites,initPrebas$maxThin,9))
initPrebas$thinning[,,9] <- -999
initPrebas$thinning[1:pine$nSites,,1:8] <- pine$thinning
initPrebas$thinning[(pine$nSites+1):(pine$nSites+spruce$nSites),1:spruce$maxThin,1:8] <- spruce$thinning

old$siteInfo[,2]<- old$siteInfo[,2] + pine$nClimID + spruce$nClimID+nfi$nClimID
nfi$siteInfo[,2] <- nfi$siteInfo[,2] + pine$nClimID + spruce$nClimID
spruce$siteInfo[,2] <- spruce$siteInfo[,2] + pine$nClimID
old$siteInfo[,1]<- old$siteInfo[,1] + pine$nSites + spruce$nSites+nfi$nSites
nfi$siteInfo[,1] <- nfi$siteInfo[,1] + pine$nSites + spruce$nSites
spruce$siteInfo[,1] <- spruce$siteInfo[,1] + pine$nSites

initPrebas$siteInfo <- rbind(pine$siteInfo,spruce$siteInfo,nfi$siteInfo,old$siteInfo[,1:9])   

initPrebas$nThinning <- c(pine$nThinning,spruce$nThinning,nfi$nThinning,old$nThinning)

initPrebas$fAPAR <- matrix(0.7,initPrebas$nSites,initPrebas$maxYears)
old$initClearcut[,5]<-10
initPrebas$initClearcut <- rbind(pine$initClearcut,spruce$initClearcut,nfi$initClearcut,old$initClearcut)
initPrebas$fixBAinitClarcut <- rep(1,initPrebas$nSites)
initPrebas$initCLcutRatio <- matrix(0,initPrebas$nSites,3) #what's the 3 here, nlayers

initPrebas$ETSy <- matrix(0,initPrebas$nClimID,initPrebas$maxYears)
initPrebas$ETSy[1:pine$nClimID,] <- pine$ETSy
initPrebas$ETSy[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears] <- spruce$ETSy
initPrebas$ETSy[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears] <- nfi$ETSy
initPrebas$ETSy[(1+pine$nClimID+spruce$nClimID+nfi$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID+old$nClimID),1:old$maxYears] <- old$ETSy

initPrebas$P0y <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,2))
initPrebas$P0y[1:pine$nClimID,,] <- pine$P0y
initPrebas$P0y[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,] <- spruce$P0y
initPrebas$P0y[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,] <- nfi$P0y
initPrebas$P0y[(1+pine$nClimID+spruce$nClimID+nfi$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID+old$nClimID),1:old$maxYears,] <- old$P0y

initPrebas$multiInitVar <- array(0,dim=c(initPrebas$nSites,7,21)) #here the 3 is layer
initPrebas$multiInitVar[1:pine$nSites,,1] <- pine$multiInitVar
initPrebas$multiInitVar[(pine$nSites+1):(pine$nSites + spruce$nSites),,1] <- spruce$multiInitVar
initPrebas$multiInitVar[(pine$nSites+1+ spruce$nSites):(nfi$nSites + pine$nSites + spruce$nSites),,] <- nfi$multiInitVar
initPrebas$multiInitVar[(pine$nSites+1+ spruce$nSites+nfi$nSites):(old$nSites+nfi$nSites + pine$nSites + spruce$nSites),,] <- old$multiInitVar

initPrebas$weather <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,365,5)) #5: number of variables 
initPrebas$weather[1:pine$nClimID,,,] <- pine$weather
initPrebas$weather[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,,] <- spruce$weather
initPrebas$weather[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,,] <- nfi$weather
initPrebas$weather[(1+pine$nClimID+spruce$nClimID+nfi$nClimID):(old$nClimID+nfi$nClimID+pine$nClimID+spruce$nClimID),1:old$maxYears,,] <- old$weather

initPrebas$soilC <- array(0,dim=c(initPrebas$nSites,initPrebas$maxYears,5,3,3))

initPrebas$weatherYasso <- array(0,dim=c(initPrebas$nClimID,initPrebas$maxYears,3))
# initPrebas$weatherYasso[1:pine$nClimID,,] <- pine$weatherYasso
# initPrebas$weatherYasso[(1+pine$nClimID):(pine$nClimID+spruce$nClimID),1:spruce$maxYears,] <- spruce$weatherYasso
# initPrebas$weatherYasso[(1+pine$nClimID+spruce$nClimID):(nfi$nClimID+pine$nClimID+spruce$nClimID),1:nfi$maxYears,] <- nfi$weatherYasso

initPrebas$soilCtot <- matrix(0,initPrebas$nSites,initPrebas$maxYears)

initPrebas$defaultThin <- rep(0,initPrebas$nSites)

initPrebas$ClCut <- rep(0,initPrebas$nSites)

initPrebas$inDclct <- matrix(0,initPrebas$nSites,3)
initPrebas$inAclct <- matrix(0,initPrebas$nSites,3)

initPrebas$dailyPRELES <- array(0,dim=c(initPrebas$nSites,(initPrebas$maxYears*365),3))

initPrebas$yassoRun <- rep(0,initPrebas$nSites)

initPrebas$PREBASversion <- rep(0,initPrebas$nSites)

spruce_B[,1] <- spruce_B[,1] + pine$nSites
spruce_D[,1] <- spruce_D[,1] + pine$nSites
spruce_H[,1] <- spruce_H[,1] + pine$nSites
spruce_Hc[,1] <- spruce_Hc[,1] + pine$nSites
spruce_V[,1] <- spruce_V[,1] + pine$nSites

nfi_B[,1] <- nfi_B[,1] + pine$nSites + spruce$nSites
nfi_D[,1] <- nfi_D[,1] + pine$nSites + spruce$nSites
nfi_H[,1] <- nfi_H[,1] + pine$nSites + spruce$nSites

old_B[,1] <- old_B[,1] + pine$nSites + spruce$nSites +nfi$nSites
old_D[,1] <- old_D[,1] + pine$nSites + spruce$nSites +nfi$nSites
old_H[,1] <- old_H[,1] + pine$nSites + spruce$nSites +nfi$nSites
old_Hc[,1] <- old_Hc[,1] + spruce$nSites + pine$nSites +nfi$nSites
old_V[,1] <- old_V[,1] + spruce$nSites + pine$nSites +nfi$nSites

outdata_B <- rbind(pine_B,spruce_B,nfi_B,old_B)
outdata_D <- rbind(pine_D,spruce_D,nfi_D,old_D)
outdata_H <- rbind(pine_H,spruce_H,nfi_H,old_H)
outdata_Hc <- rbind(spruce_Hc,old_Hc)  # what about the other plots? pine and nfi? 
outdata_V <- rbind(pine_V,spruce_V,old_V)

obs_B <- c(obs_pB,obs_sB,obs_nB,obs_oB)
obs_D <- c(obs_pD,obs_sD,obs_nD,obs_oD)
obs_H <- c(obs_pH,obs_sH,obs_nH,obs_oH)
obs_Hc <- c(obs_sHc,obs_oHc)
obs_V <- c(obs_pV,obs_sV,obs_oV)

CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,initPrebas$nClimID,initPrebas$maxYears*365)

for(i in 1:initPrebas$nClimID){
  PAR[i,] <- as.vector(t(initPrebas$weather[i,,,1]))
  TAir[i,] <- as.vector(t(initPrebas$weather[i,,,2]))
  Precip[i,] <- as.vector(t(initPrebas$weather[i,,,4]))
  VPD[i,] <- as.vector(t(initPrebas$weather[i,,,3]))
  CO2[i,] <- as.vector(t(initPrebas$weather[i,,,5]))
}

multiThin <- initPrebas$thinning
nThinning <- initPrebas$nThinning
initPrebas$siteInfo <- cbind(initPrebas$siteInfo,matrix(pPREL[1:3],nrow(initPrebas$siteInfo),3,byrow = T))
pCROBAS <- pCROB
#pCROBAS[12,] <- pCROBAS[12,] - 1 ##why?

initPrebas <- InitMultiSite(nYearsMS = initPrebas$nYears,
                            siteInfo=as.matrix(initPrebas$siteInfo),
                            pCROBAS = pCROBAS,
                            # litterSize = litterSize,#pAWEN = parsAWEN,
                            defaultThin=initPrebas$defaultThin,
                            ClCut = as.integer(initPrebas$ClCut),
                            multiInitVar = as.array(initPrebas$multiInitVar),
                            PAR = as.matrix(PAR),
                            TAir= as.matrix(TAir),
                            VPD= as.matrix(VPD),
                            Precip= as.matrix(Precip),
                            CO2= as.matrix(CO2),
                            yassoRun = initPrebas$yassoRun,#lukeRuns = initPrebas$lukeRuns,
                            # initCLcutRatio = initCLcutRatio
                            multiThin = as.array(multiThin),
                            multiNthin = as.matrix(nThinning)
)
save(initPrebas,obs_H,obs_D,obs_B,obs_V,obs_Hc,
     outdata_H,outdata_D,outdata_B,outdata_V,outdata_Hc,
     file="inputs/initPrebasAll.rdata")


# parSel <- c(1:18,31:34,38)
# nparCROB <- length(parSel)
# load("C:/Users/minunno/Documents/research/calibrationNew/nfi/chains/resDE_c24_n0.RData")
# psetNFI <- MAP(resDE)[[1]][1:(nparCROB*3)]
# load("C:/Users/minunno/Documents/research/calibrationNew/pine/chains/resDE_c15_n6x.RData")
# psetPine <- MAP(resDE)[[1]][1:nparCROB]
# load("C:/Users/minunno/Documents/research/calibrationNew/spruce/chains/resDE_c24_n2.RData")
# psetSpruce <- MAP(resDE)[[1]][1:nparCROB]
# 
# pars <- read.csv("input/parameters.csv")
# nPar <- nrow(pars)
# startValue <- matrix(runif(nPar*3,pars$parmin,pars$parmax),3,nPar,byrow = T)
# 
# startValue[1,1:(nparCROB*3)] <- psetNFI
# startValue[2,1:(nparCROB*3)] <- psetNFI
# startValue[2,1:(nparCROB)] <- psetPine
# startValue[2,(nparCROB+1):(nparCROB*2)] <- psetSpruce
# 
# save(startValue,file="input/startValue.rdata")
