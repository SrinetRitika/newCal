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
    TAir <- aperm(init$weather[,,,2], c(1,3,2))
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
  TAir <- aperm(init$weather[,,,2], c(1,3,2))
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
  


  
  
  
  multiSiteInit=init_set2
                      fertThin = 0
                          nYearsFert = 20
                          oldLayer=0
    ###initialize siteType
    multiSiteInit$multiOut[,,3,,1] <- array(multiSiteInit$siteInfo[,3],
                                            dim=c(multiSiteInit$nSites,
                                                  multiSiteInit$maxYears,
                                                  multiSiteInit$maxNlayers))
    
    ###calculate ETSmean based on a moving average window
    ETSx <- cbind(multiSiteInit$ETSstart,multiSiteInit$ETSy)
    ETSmean <- t(apply(ETSx,1,calETSmean))
    for(ijj in 1:multiSiteInit$maxNlayers){
      for(ijx in 1:multiSiteInit$nClimID){
        siteXs <- which(multiSiteInit$siteInfo[,2]==ijx)
        multiSiteInit$multiOut[siteXs,,5,ijj,2] <- rep(ETSmean[ijx,],each=length(siteXs))
      }
    }
    
    #initialize alfar
    if(is.null(multiSiteInit$pCN_alfar)){
      for(ijj in 1:multiSiteInit$maxNlayers){
        siteXs <- which(multiSiteInit$multiInitVar[,1,ijj] %in% 1:ncol(multiSiteInit$pCROBAS))
        multiSiteInit$multiOut[siteXs,,3,ijj,2] =
          matrix(multiSiteInit$pCROBAS[cbind((20+pmin(multiSiteInit$siteInfo[,3],5))[siteXs],
                                             multiSiteInit$multiInitVar[siteXs,1,ijj])],
                 length(siteXs),multiSiteInit$maxYears)
      }
    }else{
      multiSiteInit$pCROBAS[21:22,] <- multiSiteInit$pCN_alfar
      multiSiteInit$pCROBAS[23,] <- -999
      for(ijj in 1:multiSiteInit$maxNlayers){
        siteXs <- which(multiSiteInit$multiInitVar[,1,ijj] %in% 1:ncol(multiSiteInit$pCROBAS))
        alfar_p1 <- 
          matrix(multiSiteInit$pCN_alfar[1,multiSiteInit$multiInitVar[siteXs,1,ijj]],
                 length(siteXs),multiSiteInit$maxYears)
        alfar_p2 <- 
          matrix(multiSiteInit$pCN_alfar[2,multiSiteInit$multiInitVar[siteXs,1,ijj]],
                 length(siteXs),multiSiteInit$maxYears)
        CNratioSites <- CNratio(ETSmean[multiSiteInit$siteInfo[siteXs,2],],
                                multiSiteInit$multiOut[siteXs,,3,ijj,1]
                                ,pars=multiSiteInit$pECMmod[6:8])
        multiSiteInit$multiOut[siteXs,,3,ijj,2] <-  alfar_p1* exp(alfar_p2*CNratioSites) 
      }
    }
    
    if(oldLayer==1){
      multiSiteInit <- addOldLayer(multiSiteInit)
    }
    
    ####avoid species = 0  replace with species 1 when layer is empty
    multiSiteInit$multiInitVar[,1,][which(multiSiteInit$multiInitVar[,1,]==0)] <- 1
    multiSiteInit$multiOut[,,4,,1][which(multiSiteInit$multiOut[,,4,,1]==0)] = 1
    
    prebas <- .Fortran("multiPrebas",
                       multiOut = as.array(multiSiteInit$multiOut),
                       nSites = as.integer(multiSiteInit$nSites),
                       nClimID = as.integer(multiSiteInit$nClimID),
                       nLayers = as.integer(multiSiteInit$nLayers),######
                       maxYears = as.integer(multiSiteInit$maxYears),
                       maxThin = as.integer(multiSiteInit$maxThin),
                       nYears = as.integer(multiSiteInit$nYears),
                       thinning=as.array(multiSiteInit$thinning),
                       pCROBAS = as.matrix(multiSiteInit$pCROBAS),    ####
                       allSp = as.integer(multiSiteInit$allSp),       ####
                       siteInfo = as.matrix(multiSiteInit$siteInfo[,c(1:7,10:12)]),  ####
                       maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                       nThinning=as.integer(multiSiteInit$nThinning),
                       fAPAR=as.matrix(multiSiteInit$fAPAR),
                       initClearcut=as.matrix(multiSiteInit$initClearcut),
                       fixBAinitClearcut = as.double(multiSiteInit$fixBAinitClarcut),
                       initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                       ETSy=as.matrix(multiSiteInit$ETSy),
                       P0y=as.array(multiSiteInit$P0y),
                       multiInitVar=as.array(multiSiteInit$multiInitVar),
                       weather=as.array(multiSiteInit$weather),
                       DOY= as.integer(multiSiteInit$DOY),
                       pPRELES=as.double(multiSiteInit$pPRELES),
                       etmodel=as.integer(multiSiteInit$etmodel),
                       soilC = as.array(multiSiteInit$soilC),
                       pYASSO=as.double(multiSiteInit$pYASSO),
                       pAWEN = as.matrix(multiSiteInit$pAWEN),
                       weatherYasso = as.array(multiSiteInit$weatherYasso),
                       litterSize = as.array(multiSiteInit$litterSize),
                       soilCtot = as.matrix(multiSiteInit$soilCtot),
                       defaultThin=as.double(multiSiteInit$defaultThin),
                       ClCut=as.double(multiSiteInit$ClCut),
                       energyCut=as.double(multiSiteInit$energyCut),
                       inDclct=as.matrix(multiSiteInit$inDclct),
                       inAclct=as.matrix(multiSiteInit$inAclct),
                       dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                       yassoRun=as.double(multiSiteInit$yassoRun),
                       multiEnergyWood = as.array(multiSiteInit$multiEnergyWood),
                       tapioPars = as.array(multiSiteInit$tapioPars),
                       thdPer=as.double(multiSiteInit$thdPer),
                       limPer=as.double(multiSiteInit$limPer),
                       ftTapioPar = as.array(multiSiteInit$ftTapioPar),
                       tTapioPar = as.array(multiSiteInit$tTapioPar),
                       GVout = as.array(multiSiteInit$GVout),
                       GVrun = as.integer(multiSiteInit$GVrun),
                       thinInt=as.double(multiSiteInit$thinInt),
                       fertThin = as.integer(fertThin),
                       flagFert = as.integer(0),
                       nYearsFert = as.integer(nYearsFert),
                       oldLayer=as.integer(oldLayer),
                       mortMod=as.double(multiSiteInit$mortMod),
                       ECMmod=as.integer(multiSiteInit$ECMmod),
                       pECMmod=as.double(multiSiteInit$pECMmod),
                       ETSstart=as.double(multiSiteInit$ETSstart)
    )
    dimnames(prebas$multiOut) <- dimnames(multiSiteInit$multiOut)
    dimnames(prebas$multiInitVar) <- dimnames(multiSiteInit$multiInitVar)
    names(prebas$siteInfo) <- names(multiSiteInit$siteInfo)
    
    class(prebas) <- "multiPrebas"
    return(prebas)
  }
  
  