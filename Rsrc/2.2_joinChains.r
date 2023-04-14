.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library(runjags)
library(coda)

setwd("/scratch/project_2000994/calibrations/srinet/newCal/")

load("chains/calOut_0.1.rdata")
calSets <- 0:3
###settings
npar <- calOut$setup$numPars
indRun <-20 #number of independent calibration runs
thin=100
# lChain <- dim(calOut$chain[[1]])[1]
# seqX <- seq(thin,lChain,by=thin)

pChain <- mcmc.list()
pMAP <- NULL
for(i in 1:indRun){
  for(ij in calSets){
    load(paste0("chains/calOut_",ij,".",i,".rdata"))
    lChain <- dim(calOut$chain[[1]])[1]
    seqX <- seq(thin,lChain,by=thin)
    if(ij==calSets[1]){
      chain1 <- calOut$chain[[1]][seqX,]
      chain2 <- calOut$chain[[2]][seqX,]
      chain3 <- calOut$chain[[3]][seqX,]
    }else{
      chain1 <- rbind(chain1,calOut$chain[[1]][seqX,])
      chain2 <- rbind(chain2,calOut$chain[[2]][seqX,])
      chain3 <- rbind(chain3,calOut$chain[[3]][seqX,])
    }
  }
  #create a list of 3 chains. Note those 3 chains are the 3 parallel chains that are used in the DE (differential evolution algorithms)
  mcmcList <- mcmc.list(mcmc(chain1),mcmc(chain2),mcmc(chain3))
  
  #combine the 3 parallel chain in one unique chain
  pChain[[i]] <- combine.mcmc(mcmc.objects=mcmcList)
  
  LPind <- dim(pChain[[i]])[2]-2 #index to identify the prior in the chains
  pMAPx <- pChain[[i]][which.max(pChain[[i]][,LPind]),] ###check which one is the maximum a posteriori parameter vector
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

####remove LL,LP, LPrior from the chains
pChain2 <- list()
for(i in 1:indRun) pChain2[i] <- pChain[i][,1:npar]

### combine the independent calibrations
# Note that you first you need to identify and remove those chains that got stuck in a local maxima or clearly did not converged.
# to do that look at the marginal trace plots "outCal/tracePlots.pdf" starting from the loglikelihhod

# combined chains
pChainComb <- list()
xx <- mcmc.list(pChain2[1:5])
pChainComb[[1]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[6:10])
pChainComb[[2]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[11:15])
pChainComb[[3]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[16:20])
pChainComb[[4]] <- combine.mcmc(xx)

gelman.diag(pChainComb,multivariate = T)


###filter the chains
###find out the chains that got stuck in local maxima looking at the loglikelihood
MAPx <- rep(0,indRun)
for(i in 1:indRun) MAPx[i] <- max(pChain[[i]][which(is.na(pChain[[i]][,(npar+1)])==F),(npar+1)])

indX <- sort.int(MAPx, decreasing = T, index.return = TRUE)
filtChain <- indX$ix[-(1:6)]

# combined chains
pChainCombFilt <- list()
set1 <- seq(1,length(filtChain),by=2)
set2 <- seq(2,length(filtChain),by=2)


xx <- mcmc.list(pChain2[filtChain[set1]])
pChainCombFilt[[1]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[filtChain[set2]])
pChainCombFilt[[2]] <- combine.mcmc(xx)

gelman.diag(pChainCombFilt,multivariate = T)


tracePlot(pChainCombFilt)


