library(Rprebasso)

for(i in 1:4){

    initName <- paste0("init_set",i)
    
    load(paste0("inputs/",initName,".rdata"))
    
    toMem <- setdiff(ls(),c("initName","i"))
    
    mortMod <- 3
    GVrun <- 1
    ECMmod=1
    
    init <- get(initName)
    # init_set1$
    dims <- c(dim(init$weather)[1],dim(init$weather)[2]*dim(init$weather)[3])
    PAR <- aperm(init$weather[,,,1], c(1,3,2))
    dim(PAR) <- dims
    TAir <- aperm(init$weather[,,,4], c(1,3,2))
    dim(TAir) <- dims
    VPD <- aperm(init$weather[,,,3], c(1,3,2))
    dim(VPD) <- dims
    Precip <- aperm(init$weather[,,,4], c(1,3,2))
    dim(Precip) <- dims
    CO2 <- aperm(init$weather[,,,5], c(1,3,2))
    dim(CO2) <- dims
    
    
    
    
    # 
    # siteX=23
    # yearX=21
    # plot(PAR[siteX,((yearX-1)*365+1):(yearX*365)])
    # points(init$weather[siteX,yearX,,1],col=2,pch=20)
    
    init <- InitMultiSite(nYearsMS = init$nYears,
                             pCROBAS = pCROB,
                             pHcMod = pHcM,
                             pPRELES = pPREL,
                             pYASSO = pYAS,
                             pAWEN = parsAWEN,
                             siteInfo = init$siteInfo,
                             multiInitVar = init$multiInitVar,
                             multiThin = init$thinning,
                             multiNthin = init$nThinning,
                             multiInitClearCut = init$initClearcut,
                             initCLcutRatio = init$initCLcutRatio,
                             PAR = PAR,
                             TAir = TAir,
                             VPD = VPD,
                             Precip = Precip,
                             CO2 = CO2,
                             defaultThin = 0,
                             ClCut = 0,
                             yassoRun = 0,
                             GVrun = GVrun,
                             mortMod = mortMod,
                             ECMmod = ECMmod,
                             pCN_alfar = parsCN_alfar)
    
    assign(initName,init)
    save(list=toMem, file=paste0("inputs/",initName,"_newVersion.rdata"))
    
    print(initName)    
    rm(list=ls())
    gc()

}


library(Rprebasso)


  initName <- "init_set5Flux"
  
  load(paste0("inputs/",initName,".RData"))
  
  toMem <- setdiff(ls(),c("initName","i"))
  
  mortMod <- 3
  GVrun <- 1
  ECMmod=1
  
  init <- get(initName)
  # init_set1$
  dims <- c(dim(init$weather)[1],dim(init$weather)[2]*dim(init$weather)[3])
  PAR <- aperm(init$weather[,,,1], c(1,3,2))
  dim(PAR) <- dims
  TAir <- aperm(init$weather[,,,4], c(1,3,2))
  dim(TAir) <- dims
  VPD <- aperm(init$weather[,,,3], c(1,3,2))
  dim(VPD) <- dims
  Precip <- aperm(init$weather[,,,4], c(1,3,2))
  dim(Precip) <- dims
  CO2 <- aperm(init$weather[,,,5], c(1,3,2))
  dim(CO2) <- dims
  
  
  
  
  # 
  # siteX=23
  # yearX=21
  # plot(PAR[siteX,((yearX-1)*365+1):(yearX*365)])
  # points(init$weather[siteX,yearX,,1],col=2,pch=20)
  
  init <- InitMultiSite(nYearsMS = init$nYears,
                        pCROBAS = pCROB,
                        pHcMod = pHcM,
                        pPRELES = pPREL,
                        pYASSO = pYAS,
                        pAWEN = parsAWEN,
                        siteInfo = init$siteInfo,
                        multiInitVar = init$multiInitVar,
                        multiThin = init$thinning,
                        multiNthin = init$nThinning,
                        multiInitClearCut = init$initClearcut,
                        initCLcutRatio = init$initCLcutRatio,
                        PAR = PAR,
                        TAir = TAir,
                        VPD = VPD,
                        Precip = Precip,
                        CO2 = CO2,
                        defaultThin = 0,
                        ClCut = 0,
                        yassoRun = 0,
                        GVrun = GVrun,
                        mortMod = mortMod,
                        ECMmod = ECMmod,
                        pCN_alfar = parsCN_alfar)
  
  assign(initName,init)
  save(list=toMem, file=paste0("inputs/",initName,"_newVersion.rdata"))
  
  print(initName)    
  rm(list=ls())
  gc()
  

