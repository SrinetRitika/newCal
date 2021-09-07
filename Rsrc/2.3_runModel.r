vLocal <- FALSE
if(!vLocal){
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
library(ggplot2)
library(Metrics)
library(data.table)
# library(Rprebasso)
library(BayesianTools)
library(ggpubr)

devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

# setwd("C:/Users/minunno/Documents/research/calibrationNew/all")
setwd("/scratch/project_2000994/calibrations/all")


load("inputs/initPrebas.rdata")

load("outCal/pMAP.rdata")
# fileNamelastCal <- paste0("chains/calOut_n8.RData")
# # load("input/startValue.rdata")
# if(!exists("pMAP")){
#   load(fileNamelastCal)
#   pMAP <- MAP(calOut)[[1]]
# }

parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)

initPrebas$pCROBAS <- pCROB
initPrebas$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
initPrebas$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
initPrebas$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

# pCROBAS <- initPrebas$pCROBAS 
# save(pCROBAS,file=paste0(unlist(strsplit(fileNamelastCal,"[.]"))[1],"_pMAP.rdata"))

startX <- Sys.time()
modOut <- multiPrebas(initPrebas)
endX <- Sys.time()
timeX = endX- startX
print(timeX)

# par(mfrow=c(2,2))
# outdata_V2 <- outdata_V; outdata_V2[,5] <- 2
# outdata_B2 <- outdata_B; outdata_B2[,5] <- 2
# out_Vold <- modOut$multiOut[outdata_V] + modOut$multiOut[outdata_V2]
# out_Hold <- modOut$multiOut[outdata_H]
# out_Dold <- modOut$multiOut[outdata_D]
# out_Bold <- modOut$multiOut[outdata_B] + modOut$multiOut[outdata_B2]

# load("calOut_c15_n2.RData")
# sp=1
# psetNew <- MAP(calOut)[[1]]
# parSel <- c(1:18,31:34,38)
# nparCROB <- length(parSel)
# initPrebas$pCROBAS[parSel,sp] <- psetNew[1:nparCROB]
#  
# startX <- Sys.time()
# modOut <- multiPrebas(initPrebas)
# endX <- Sys.time()
# timeX = endX- startX
# print(timeX)

pgePsites <- 1:657
pgeSPsites <- 658:785
nfisites <- 786:922
outdata_V2 <- outdata_V; outdata_V2[,5] <- 2
outdata_B2 <- outdata_B; outdata_B2[,5] <- 2
out_Vnew <- modOut$multiOut[outdata_V] + modOut$multiOut[outdata_V2]
out_Hnew <- modOut$multiOut[outdata_H]
out_Hcnew <- modOut$multiOut[outdata_Hc]
out_Dnew <- modOut$multiOut[outdata_D]
out_Bnew <- modOut$multiOut[outdata_B] + modOut$multiOut[outdata_B2]

speciesNam <- c("pine","spruce","birch")
spVind = outdata_V; spVind[,3] <- 4; spVind[,5] <- 1
allData <- data.table(obs = obs_V,sim=out_Vnew,site=outdata_V[,1],species=speciesNam[modOut$multiOut[spVind]],var="V")
spBind = outdata_B; spBind[,3] <- 4; spBind[,5] <- 1
allData <- rbind(allData,data.table(obs = obs_B,sim=out_Bnew,site=outdata_B[,1],species=speciesNam[modOut$multiOut[spBind]],var="B"))
spHind = outdata_H; spHind[,3] <- 4; spHind[,5] <- 1
allData <- rbind(allData,data.table(obs = obs_H,sim=out_Hnew,site=outdata_H[,1],species=speciesNam[modOut$multiOut[spHind]],var="H"))
spHcind = outdata_Hc; spHcind[,3] <- 4; spHcind[,5] <- 1
allData <- rbind(allData,data.table(obs = obs_Hc,sim=out_Hcnew,site=outdata_Hc[,1],species=speciesNam[modOut$multiOut[spHcind]],var="Hc"))
spDind = outdata_D; spDind[,3] <- 4; spDind[,5] <- 1
allData <- rbind(allData,data.table(obs = obs_D,sim=out_Dnew,site=outdata_D[,1],species=speciesNam[modOut$multiOut[spDind]],var="D"))
allData[site %in% pgePsites,dataset:="PGEpine"]
allData[site %in% pgeSPsites,dataset:="PGEspruce"]
allData[site %in% nfisites,dataset:="NFI"]

pdf(file="outCal/outPlots.pdf")
pV <- ggplot(allData[var=="V"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("V")
pB <- ggplot(allData[var=="B"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("B")
pD <- ggplot(allData[var=="D"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("D")
pH <- ggplot(allData[var=="H"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("H")
pHc <- ggplot(allData[var=="Hc"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("Hc")
ggarrange(pB,pV,pD,pH,pHc,common.legend = F)
print(pV)
print(pB)
print(pD)
print(pH)
print(pHc)


###pgePIne
pV <- ggplot(allData[var=="V" & dataset=="PGEpine"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("V PGEpine")
pB <- ggplot(allData[var=="B" & dataset=="PGEpine"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("B PGEpine")
pD <- ggplot(allData[var=="D" & dataset=="PGEpine"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("D PGEpine")
pH <- ggplot(allData[var=="H" & dataset=="PGEpine"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("H PGEpine")
pHc <- ggplot(allData[var=="Hc" & dataset=="PGEpine"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("Hc PGEpine")
ggarrange(pB,pV,pD,pH,common.legend = F)
print(pV)
print(pB)
print(pD)
print(pH)

###pgeSpruce
pV <- ggplot(allData[var=="V" & dataset=="PGEspruce"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("V PGEspruce")
pB <- ggplot(allData[var=="B" & dataset=="PGEspruce"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("B PGEspruce")
pD <- ggplot(allData[var=="D" & dataset=="PGEspruce"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("D PGEspruce")
pH <- ggplot(allData[var=="H" & dataset=="PGEspruce"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("H PGEspruce")
pHc <- ggplot(allData[var=="Hc" & dataset=="PGEspruce"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("Hc PGEspruce")
ggarrange(pB,pV,pD,pH,pHc,common.legend = F)
print(pV)
print(pB)
print(pD)
print(pH)
print(pHc)

###NFI
pB <- ggplot(allData[var=="B" & dataset=="NFI"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("B NFI")
pD <- ggplot(allData[var=="D" & dataset=="NFI"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("D NFI")
pH <- ggplot(allData[var=="H" & dataset=="NFI"],aes(x=obs,y=sim,col=species)) + 
  geom_point() + geom_abline() + ggtitle("H NFI")
ggarrange(pB,pD,pH,common.legend = F)
print(pB)
print(pD)
print(pH)
dev.off()


pdf("pDist.pdf",30,30)
 marginalPlot(calOut,singlePanel = F,prior = F)#,xrange = t(xrange))
 correlationPlot(calOut,singlePanel = F,scaleCorText = FALSE, thin = 100, start = 20000)
dev.off()

if(T){
  pdf(file="outPlots2.pdf")
  findDataSt <- function(st,initPrebas,outdata,obs,out,varX){
    stX <- which(initPrebas$siteInfo[,3]==st)
    dataX <- which(outdata[,1] %in% stX)
    obs[dataX]
    out[dataX]
    par(mfrow=c(1,2))
    plot(obs,out,main=paste(varX,"st =",st),ylab="simulated",xlab="observed")
    points(obs[dataX],out[dataX],col=2,pch=20)
    abline(0,1)
    rmseX <- rmse(obs,out)
    legend("bottomright", paste0("rmse ",round(rmseX,2)),bty="n") 
    residX <- obs[dataX]-out[dataX]
    hist(residX,main=NULL,xlab="residuals")
    metricX <- c(min(residX),mean(residX),median(residX),max(residX))
    legend("topleft",legend=round(metricX,1),bty="n")
  }
  
  findDataSt(2,initPrebas,outdata_H,obs_H,out_Hnew,"H")
  findDataSt(3,initPrebas,outdata_H,obs_H,out_Hnew,"H")
  findDataSt(4,initPrebas,outdata_H,obs_H,out_Hnew,"H")
  findDataSt(5,initPrebas,outdata_H,obs_H,out_Hnew,"H")
  
  findDataSt(2,initPrebas,outdata_D,obs_D,out_Dnew,"D")
  findDataSt(3,initPrebas,outdata_D,obs_D,out_Dnew,"D")
  findDataSt(4,initPrebas,outdata_D,obs_D,out_Dnew,"D")
  findDataSt(5,initPrebas,outdata_D,obs_D,out_Dnew,"D")
  
  findDataSt(2,initPrebas,outdata_B,obs_B,out_Bnew,"B")
  findDataSt(3,initPrebas,outdata_B,obs_B,out_Bnew,"B")
  findDataSt(4,initPrebas,outdata_B,obs_B,out_Bnew,"B")
  findDataSt(5,initPrebas,outdata_B,obs_B,out_Bnew,"B")
  
  findDataSt(2,initPrebas,outdata_V,obs_V,out_Vnew,"V")
  findDataSt(3,initPrebas,outdata_V,obs_V,out_Vnew,"V")
  findDataSt(4,initPrebas,outdata_V,obs_V,out_Vnew,"V")
  findDataSt(5,initPrebas,outdata_V,obs_V,out_Vnew,"V")
  
  
  pX <- function(out,varX,siteX,spX,obs,outdata,outX){
    startYear <- modOut$multiOut[siteX,1,7,1,1]
    nYears = out$nYears[siteX]
    dataX <- data.table(out$multiOut[siteX,1:nYears,varX,spX,1])
    indX <- which(outdata[,1]==siteX & outdata[,3]==varX)
    if(length(indX)>0){
      yearData <- outdata[indX,2]
      sim <- outX[indX]
      ggp <- ggplot() + 
        geom_line(data = dataX, aes(x=((1:nYears)+startYear), y=V1,col="sim")) +
        geom_point(aes(x=(yearData+startYear), y=sim,col="sim")) +
        geom_point(aes(x=(yearData+startYear), y=obs[indX],col="obs")) + 
        ggtitle(paste(varNames[varX],"site ",siteX))+
        ylab("") + xlab("age") +
        theme(legend.title = element_blank())
    } else{
      ggp <- ggplot() + 
        geom_line(data = dataX, aes(x=((1:nYears)+startYear), y=V1,col="sim")) +
        ggtitle(paste(varNames[varX],"site ",siteX)) +
        ylab("") + xlab("age") +
        theme(legend.title = element_blank())
    }
    
  }
  
  nSites <- modOut$nSites
  
  for(siteX in 1:nSites){
    p1 <- pX(modOut,11,siteX,1,obs_H,outdata_H,out_Hnew)
    p2 <- pX(modOut,12,siteX,1,obs_D,outdata_D,out_Dnew)
    p3 <- pX(modOut,13,siteX,1,obs_B,outdata_B,out_Bnew)
    p4 <- pX(modOut,14,siteX,1,NULL,NULL,NULL)
    p5 <- pX(modOut,30,siteX,1,obs_V,outdata_V,out_Vnew)
    p6 <- pX(modOut,17,siteX,1,NULL,NULL,NULL)
    
    print(ggarrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3))
  }  
  dev.off()
}


