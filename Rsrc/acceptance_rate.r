.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library(runjags)
library(coda)

setwd("/scratch/project_2000994/calibrations/srinet/newCal/")

load("chains/calOut_29.1.rdata")
# parmod=MAP(calOut)$parametersMAP
# pdf("outCal/marginalPlots_25_5.pdf")
# tracePlot(calOut)
# marginalPlot(calOut,prior=T)
# gelmanDiagnostics(calOut)
# dev.off()

calSet <- 29
###settings
npar <- calOut$setup$numPars
indRun <-1 #number of independent calibration runs
arr_lp <- array(NA, dim=c(indRun,4,3))
dimnames(arr_lp) <- list(NULL, c("max","min","range","acc. rate"))
pChain <- mcmc.list()

for(i in 1:indRun){
    load(paste0("chains/calOut_",calSet,".",i,".rdata"))
    for(k in 1:3){
      arr_lp[i,1,k] <- max(calOut$chain[[k]][,(npar+1)])
      arr_lp[i,2,k] <- min(calOut$chain[[k]][,(npar+1)])
      arr_lp[i,3,k] <- arr_lp[i,1,k]-arr_lp[i,2,k]
      arr_lp[i,4,k] <- length(unique(calOut$chain[[k]][,npar+1]))/10000
    }
  
  lChain <- dim(calOut$chain[[1]])[1]
    chain1 <- calOut$chain[[1]]
    chain2 <- calOut$chain[[2]]
    chain3 <- calOut$chain[[3]]
 
mcmcList <- mcmc.list(mcmc(chain1),mcmc(chain2),mcmc(chain3))
  
  #combine the 3 parallel chain in one unique chain
  pChain[[i]] <- combine.mcmc(mcmc.objects=mcmcList)
}

pdf("outCal/tracePlots_29.pdf")
tracePlot(pChain)
dev.off()
