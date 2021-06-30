.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library(runjags)
library(coda)

setwd("/scratch/project_2000994/calibrations/all")

load("chains/calOut_2.1.rdata")

###settings
thin=1000
lChain <- dim(calOut$chain[[1]])[1]
seqX <- seq(thin,lChain,by=thin)


pChain <- mcmc.list()
pMAP <- NULL
for(i in 1:20){
  load(paste0("chains/calOut_1.",i,".rdata"))
  chain1 <- calOut$chain[[1]][seqX,]
  chain2 <- calOut$chain[[2]][seqX,]
  chain3 <- calOut$chain[[3]][seqX,]
  load(paste0("chains/calOut_2.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_3.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_4.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_5.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_6.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_7.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_8.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  load(paste0("chains/calOut_9.",i,".rdata"))
  chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  # load(paste0("chains/calOut_10.",i,".rdata"))
  # chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
  # chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
  # chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
  
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

