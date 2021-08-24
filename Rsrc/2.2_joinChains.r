.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library(runjags)
library(coda)

setwd("/scratch/project_2000994/calibrations/newCal")

load("chains/calOut_0.1.rdata")

###settings
thin=1000
lChain <- dim(calOut$chain[[1]])[1]
seqX <- seq(thin,lChain,by=thin)

pChain <- mcmc.list()
pMAP <- NULL
for(i in 1:20){
  for(ij in calSets){
    load(paste0("chains/calOut_",ij,".",i,".rdata"))
    chain1 <- calOut$chain[[1]][seqX,]
    chain2 <- calOut$chain[[2]][seqX,]
    chain3 <- calOut$chain[[3]][seqX,]
  }
  mcmcList <- mcmc.list(mcmc(chain1),mcmc(chain2),mcmc(chain3))
  
  pChain[[i]] <- combine.mcmc(mcmc.objects=mcmcList)

  LPind <- dim(pChain[[i]])[2]-2
  pMAPx <- pChain[[i]][which.max(pChain[[i]][,LPind]),]
  if(i ==1){
    pMAP <- pMAPx
  }else{
    if(pMAPx[LPind] > pMAP[LPind]) pMAP <- pMAPx  
  }
  
}

save(pChain, file="outCal/allChain.rdata")
save(pMAP,file="outCal/pMAP.rdata")

pdf("outCal/tracePlots.pdf")
 tracePlot(pChain)
dev.off()

