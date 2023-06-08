library(dplyr)
library(data.table)
library("readxl")
library(Rprebasso)
library(stringr)
library(tidyverse)

dataSp <- data.table(read_excel("inputs/Spruce-experiments-area_31May.xlsx", sheet = "Kokeet-5-2022"))
dataPi <- data.table(read_excel("inputs/Pine-experiments-area_31May.xlsx", sheet = "Kokeet-5-2022"))
load('C:/Users/srinetri/research/PREBAS_calibration/newCal/inputs/init_set1.rdata')
load('C:/Users/srinetri/research/PREBAS_calibration/newCal/inputs/init_set2.rdata')
load('C:/Users/srinetri/research/PREBAS_calibration/newCal/inputs/init_set3.rdata')
load('C:/Users/srinetri/research/PREBAS_calibration/newCal/inputs/init_set4.rdata')

### Spruce - init_set3
dataSp<- dataSp %>%
  mutate(index = match(str_c(ExpID,Plot), unique(str_c(ExpID,Plot)))) %>%
  arrange(index)
nSites<-length(unique(dataSp$index))
dataSp$BasalAreaKL<-NA
dataSp$NStemsKL<-NA
dataSp$VolumeKL<-NA
dataSp$quadDKL<-NA
for (i in 1:nSites){
  for (j in 1:length(unique(dataSp$Age[which(dataSp$index==i)]))){
    ageX=unique(dataSp$Age[which(dataSp$index==i)])[j]
    if (dataSp$index==i && dataSp$Age == ageX){
      dataSp$BasalAreaKL[dataSp$Segment == 'K']<-dataSp$BasalArea[dataSp$Segment == 'K']-dataSp$BasalArea[dataSp$Segment == 'LP']
      dataSp$NStemsKL[dataSp$Segment == 'K']<-dataSp$NStems[dataSp$Segment == 'K']-dataSp$NStems[dataSp$Segment == 'LP']
      dataSp$VolumeKL[dataSp$Segment == 'K']<-dataSp$Volume[dataSp$Segment == 'K']-dataSp$Volume[dataSp$Segment == 'LP']
      dataSp$quadDKL[dataSp$Segment == 'K']<- sqrt(dataSp$BasalAreaKL[dataSp$Segment == 'K']/dataSp$NStemsKL[dataSp$Segment == 'K']/pi*10000)*2
      dataSp$quadDKL[dataSp$Segment == 'J']<- sqrt(dataSp$BasalArea[dataSp$Segment == 'J']/dataSp$NStems[dataSp$Segment == 'J']/pi*10000)*2
    }
  }
}
dataSp[,year1:=min(Year),by=index]
dataSp[,simYear:=Year- year1]
dataSp[,status:=1]
dataSp[,status:=ifelse(BasalArea[Segment=="L"]==0,1,2),by=.(index, MeasurementTime)]

initdataSp<-dataSp[simYear==0 & Segment=="J"]
spX2 <- 2
baX <- initdataSp$BasalArea
dwX <- as.numeric(initdataSp$MeanD_Baweighted)
dQ <- initdataSp$quadDKL

speciesFilter <- which(init_set3$multiInitVar[,1,1]==spX2)
init<-as.data.frame(init_set3$multiInitVar[speciesFilter,,1])
init<-init %>%
  mutate(matches = pmap(list(D, BA),
                        function(dw, ba) filter(initdataSp,
                                                MeanD_Baweighted == dw,
                                                BasalArea == ba)$index),)
for(i in 1:length(init$matches))if (length(init$matches[[i]]) == 0){
  return(init$matches[[i]]<-NA)
return(init$matches[[i]])
}
init$matches<-unlist(init$matches)
indN<-which(is.na(init$matches)==F)
init_set3Key<-data.frame(Index=init$matches[indN],siteX=speciesFilter[indN])

nLayers<-max(init_set3$nLayers)
nYears<-init_set3$nYears[speciesFilter[indN]]
nSites<-length(speciesFilter[indN])

multiInitVarNew <- array(NA, dim=c(length(speciesFilter[indN]), 7,nLayers),dimnames =
                           dimnames(init_set3$multiInitVar))
siteInfoNew <- matrix(0, length(speciesFilter[indN]), 12, byrow = T,
                      dimnames = dimnames(init_set3$siteInfo))
thinningNew <- array(NA,dim=c(length(speciesFilter[indN]),4,10))

nThinningNew<-c(length(speciesFilter[indN]))

for(i in 1:length(speciesFilter[indN])){
  init_set3$multiInitVar[speciesFilter[indN][i],4,1]<-dQ[init$matches[indN][i]]
  multiInitVarNew[i,,]<-init_set3$multiInitVar[speciesFilter[indN][i],,]
  siteInfoNew[i,]<-init_set3$siteInfo[speciesFilter[indN][i],]
  thinningNew[i,,]<-init_set3$thinning[speciesFilter[indN][i],,]
  nThinningNew[i]<-init_set3$nThinning[speciesFilter[indN][i]]
}
weatherNew<-init_set3$weather[unique(siteInfoNew[,2]),,,]

climIDx <- sort(unique(siteInfoNew[,2]))
siteInfoNew[,2] <- match(siteInfoNew[,2],climIDx)
maxYears <- max(init_set3$nYears[speciesFilter[indN]])
CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)

for(i in 1:length(climIDx)){
  ij <- climIDx[i]
  PAR[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,1]))
  TAir[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,2]))
  Precip[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,4]))
  VPD[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,3]))
  CO2[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,5]))
}

obsDout <- dataSp[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=12,layer=1,status=status,observed=quadDKL)]
obsDout_s3Sp<-obsDout[is.element(obsDout$siteID,intersect(obsDout$siteID,init_set3Key$Index)),]
obsDout_s3Sp<-obsDout_s3Sp %>% mutate(siteID = coalesce(
  init_set3Key$siteX[match(siteID, init_set3Key$Index)]
))
obsBout <- dataSp[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=13,layer=1,status=status,observed=BasalAreaKL)]
obsBout_s3Sp<-obsBout[is.element(obsBout$siteID,intersect(obsBout$siteID,init_set3Key$Index)),]
obsBout_s3Sp<-obsBout_s3Sp %>% mutate(siteID = coalesce(
  init_set3Key$siteX[match(siteID, init_set3Key$Index)]
))
obsHout <- dataSp[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=11,layer=1,status=status,observed=Hdom)]
obsHout_s3Sp<-obsHout[is.element(obsHout$siteID,intersect(obsHout$siteID,init_set3Key$Index)),]
obsHout_s3Sp<-obsHout_s3Sp %>% mutate(siteID = coalesce(
  init_set3Key$siteX[match(siteID, init_set3Key$Index)]
))
obsVout <- dataSp[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=30,layer=1,status=status,observed=VolumeKL)]
obsVout_s3Sp<-obsVout[is.element(obsVout$siteID,intersect(obsVout$siteID,init_set3Key$Index)),]
obsVout_s3Sp<-obsVout_s3Sp %>% mutate(siteID = coalesce(
  init_set3Key$siteX[match(siteID, init_set3Key$Index)]
))
obsNout <- dataSp[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=17,layer=1,status=status,observed=NStemsKL)]
obsNout_s3Sp<-obsNout[is.element(obsNout$siteID,intersect(obsNout$siteID,init_set3Key$Index)),]
obsNout_s3Sp<-obsNout_s3Sp %>% mutate(siteID = coalesce(
  init_set3Key$siteX[match(siteID, init_set3Key$Index)]
))

init_set3NewSpruce<-InitMultiSite(nYearsMS = nYears,
              siteInfo = siteInfoNew,
              multiInitVar = multiInitVarNew,
              multiThin = thinningNew,
              multiNthin = nThinningNew,
              PAR = PAR,
              TAir = TAir,
              VPD = VPD,
              Precip = Precip,
              CO2 = CO2,
              defaultThin = 0,
              ClCut = 0)

init_set3NewSpruce<-save(init_set3NewSpruce,obsDout_s3Sp,obsBout_s3Sp,obsHout_s3Sp,obsVout_s3Sp,obsNout_s3Sp, file="inputs/init_set3NewSpruce.rdata")

### Pine
dataPi<- dataPi %>%
  mutate(index = match(str_c(ExpID,Plot), unique(str_c(ExpID,Plot)))) %>%
  arrange(index)
nSites<-length(unique(dataPi$index))
dataPi$BasalAreaKL<-NA
dataPi$NStemsKL<-NA
dataPi$VolumeKL<-NA
dataPi$quadDKL<-NA
for (i in 1:nSites){
  for (j in 1:length(unique(dataPi$Age[which(dataPi$index==i)]))){
    ageX=unique(dataPi$Age[which(dataPi$index==i)])[j]
    if (dataPi$index==i && dataPi$Age == ageX){
      dataPi$BasalAreaKL[dataPi$Segment == 'K']<-dataPi$BasalArea[dataPi$Segment == 'K']-dataPi$BasalArea[dataPi$Segment == 'LP']
      dataPi$NStemsKL[dataPi$Segment == 'K']<-dataPi$NStems[dataPi$Segment == 'K']-dataPi$NStems[dataPi$Segment == 'LP']
      dataPi$VolumeKL[dataPi$Segment == 'K']<-dataPi$Volume[dataPi$Segment == 'K']-dataPi$Volume[dataPi$Segment == 'LP']
      dataPi$quadDKL[dataPi$Segment == 'K']<- sqrt(dataPi$BasalAreaKL[dataPi$Segment == 'K']/dataPi$NStemsKL[dataPi$Segment == 'K']/pi*10000)*2
      dataPi$quadDKL[dataPi$Segment == 'J']<- sqrt(dataPi$BasalArea[dataPi$Segment == 'J']/dataPi$NStems[dataPi$Segment == 'J']/pi*10000)*2
    }
  }
}

dataPi[,year1:=min(Year),by=index]
dataPi[,simYear:=Year- year1]
dataPi[,status:=1]
dataPi[,status:=ifelse(BasalArea[Segment=="L"]==0,1,2),by=.(index, MeasurementTime)]


initdataPi<-dataPi[simYear==0 & Segment=="J"]
spX1 <- 1
baX <- initdataPi$BasalArea
dwX <- as.numeric(initdataPi$MeanD_Baweighted)
dQ <- initdataPi$quadDKL

speciesFilter <- which(init_set1$multiInitVar[,1,1]==spX1)
initPi1<-as.data.frame(init_set1$multiInitVar[speciesFilter,,1])
initPi1<-initPi1 %>%
  mutate(matches = pmap(list(D, BA),
                        function(dw, ba) filter(initdataPi,
                                                MeanD_Baweighted == dw,
                                                BasalArea == ba)$index),)
for(i in 1:length(initPi1$matches)) if (length(initPi1$matches[[i]]) == 0){
  initPi1$matches[[i]]<-NA
}
initPi1$matches[[196]]<-46
initPi1$matches[[197]]<-47
initPi1$matches[[198]]<-48
initPi1$matches[[199]]<-49
initPi1$matches[[200]]<-50
initPi1$matches[[201]]<-51
initPi1$matches[[262]]<-219
initPi1$matches[[263]]<-220
initPi1$matches[[264]]<-221

initPi1$matches<-unlist(initPi1$matches)
indN<-which(is.na(initPi1$matches)==F)
init_set1Key<-data.frame(Index=initPi1$matches[indN],siteX=speciesFilter[indN])

nLayers<-max(init_set1$nLayers)
nYears<-init_set1$nYears[speciesFilter[indN]]
nSites<-length(speciesFilter[indN])

multiInitVarNew <- array(NA, dim=c(length(speciesFilter[indN]), 7,nLayers),dimnames =
                           dimnames(init_set1$multiInitVar))
siteInfoNew <- matrix(0, length(speciesFilter[indN]), 12, byrow = T,
                      dimnames = dimnames(init_set1$siteInfo))
thinningNew <- array(NA,dim=c(length(speciesFilter[indN]),7,10))

nThinningNew<-c(length(speciesFilter[indN]))

for(i in 1:length(speciesFilter[indN])){
  init_set1$multiInitVar[speciesFilter[indN][i],4,1]<-dQ[initPi1$matches[indN][i]]
  multiInitVarNew[i,,]<-init_set1$multiInitVar[speciesFilter[indN][i],,]
  siteInfoNew[i,]<-init_set1$siteInfo[speciesFilter[indN][i],]
  thinningNew[i,,]<-init_set1$thinning[speciesFilter[indN][i],,]
  nThinningNew[i]<-init_set1$nThinning[speciesFilter[indN][i]]
}
weatherNew<-init_set1$weather[unique(siteInfoNew[,2]),,,]

climIDx <- sort(unique(siteInfoNew[,2]))
siteInfoNew[,2] <- match(siteInfoNew[,2],climIDx)
maxYears <- max(init_set1$nYears[speciesFilter[indN]])
CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)

for(i in 1:length(climIDx)){
  ij <- climIDx[i]
  PAR[i,] <- as.vector(t(init_set1$weather[ij,1:maxYears,,1]))
  TAir[i,] <- as.vector(t(init_set1$weather[ij,1:maxYears,,2]))
  Precip[i,] <- as.vector(t(init_set1$weather[ij,1:maxYears,,4]))
  VPD[i,] <- as.vector(t(init_set1$weather[ij,1:maxYears,,3]))
  CO2[i,] <- as.vector(t(init_set1$weather[ij,1:maxYears,,5]))
}

obsDout <- dataPi[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=12,layer=1,status=status,observed=quadDKL)]
obsDout_s1<-obsDout[is.element(obsDout$siteID,intersect(obsDout$siteID,init_set1Key$Index)),]
obsDout_s1<-obsDout_s1 %>% mutate(siteID = coalesce(
  init_set1Key$siteX[match(siteID, init_set1Key$Index)]
))
obsBout <- dataPi[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=13,layer=1,status=status,observed=BasalAreaKL)]
obsBout_s1<-obsBout[is.element(obsBout$siteID,intersect(obsBout$siteID,init_set1Key$Index)),]
obsBout_s1<-obsBout_s1 %>% mutate(siteID = coalesce(
  init_set1Key$siteX[match(siteID, init_set1Key$Index)]
))
obsHout <- dataPi[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=11,layer=1,status=status,observed=Hdom)]
obsHout_s1<-obsHout[is.element(obsHout$siteID,intersect(obsHout$siteID,init_set1Key$Index)),]
obsHout_s1<-obsHout_s1 %>% mutate(siteID = coalesce(
  init_set1Key$siteX[match(siteID, init_set1Key$Index)]
))
obsVout <- dataPi[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=30,layer=1,status=status,observed=VolumeKL)]
obsVout_s1<-obsVout[is.element(obsVout$siteID,intersect(obsVout$siteID,init_set1Key$Index)),]
obsVout_s1<-obsVout_s1 %>% mutate(siteID = coalesce(
  init_set1Key$siteX[match(siteID, init_set1Key$Index)]
))
obsNout <- dataPi[simYear>0 & Segment=="K",.(siteID=index,year=simYear,variable=17,layer=1,status=status,observed=NStemsKL)]
obsNout_s1<-obsNout[is.element(obsNout$siteID,intersect(obsNout$siteID,init_set1Key$Index)),]
obsNout_s1<-obsNout_s1 %>% mutate(siteID = coalesce(
  init_set1Key$siteX[match(siteID, init_set1Key$Index)]
))

init_set1New<-InitMultiSite(nYearsMS = nYears,
                            siteInfo = siteInfoNew,
                            multiInitVar = multiInitVarNew,
                            multiThin = thinningNew,
                            multiNthin = nThinningNew,
                            PAR = PAR,
                            TAir = TAir,
                            VPD = VPD,
                            Precip = Precip,
                            CO2 = CO2,
                            defaultThin = 0,
                            ClCut = 0)

init_set1New<-save(init_set1New,obsDout_s1,obsBout_s1,obsHout_s1,obsVout_s1,obsNout_s1, file="inputs/init_set1New.rdata")

### Pine - init_set2
speciesFilter <- which(init_set2$multiInitVar[,1,1]==spX1)
initPi2<-as.data.frame(init_set2$multiInitVar[speciesFilter,,1])
initPi2<-initPi2 %>%
  mutate(matches = pmap(list(D, BA),
                        function(dw, ba) filter(initdataPi,
                                                MeanD_Baweighted == dw,
                                                BasalArea == ba)$index),)
for(i in 1:length(initPi2$matches)) if (length(initPi2$matches[[i]]) == 0){
  initPi2$matches[[i]]<-NA
}

initPi2$matches<-unlist(initPi2$matches)
indN<-which(is.na(initPi2$matches)==F)
init_set2Key<-data.frame(Index=initPi2$matches[indN],siteX=speciesFilter[indN])

nLayers<-max(init_set2$nLayers)
nYears<-init_set2$nYears[speciesFilter[indN]]
nSites<-length(speciesFilter[indN])

multiInitVarNew <- array(NA, dim=c(length(speciesFilter[indN]), 7,nLayers),dimnames =
                           dimnames(init_set2$multiInitVar))
siteInfoNew <- matrix(0, length(speciesFilter[indN]), 12, byrow = T,
                      dimnames = dimnames(init_set2$siteInfo))
thinningNew <- array(NA,dim=c(length(speciesFilter[indN]),7,10))

nThinningNew<-c(length(speciesFilter[indN]))

for(i in 1:length(speciesFilter[indN])){
  init_set2$multiInitVar[speciesFilter[indN][i],4,1]<-dQ[initPi2$matches[indN][i]]
  multiInitVarNew[i,,]<-init_set2$multiInitVar[speciesFilter[indN][i],,]
  siteInfoNew[i,]<-init_set2$siteInfo[speciesFilter[indN][i],]
  thinningNew[i,,]<-init_set2$thinning[speciesFilter[indN][i],,]
  nThinningNew[i]<-init_set2$nThinning[speciesFilter[indN][i]]
}
weatherNew<-init_set2$weather[unique(siteInfoNew[,2]),,,]

climIDx <- sort(unique(siteInfoNew[,2]))
siteInfoNew[,2] <- match(siteInfoNew[,2],climIDx)
maxYears <- max(init_set2$nYears[speciesFilter[indN]])
CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)

for(i in 1:length(climIDx)){
  ij <- climIDx[i]
  PAR[i,] <- as.vector(t(init_set2$weather[ij,1:maxYears,,1]))
  TAir[i,] <- as.vector(t(init_set2$weather[ij,1:maxYears,,2]))
  Precip[i,] <- as.vector(t(init_set2$weather[ij,1:maxYears,,4]))
  VPD[i,] <- as.vector(t(init_set2$weather[ij,1:maxYears,,3]))
  CO2[i,] <- as.vector(t(init_set2$weather[ij,1:maxYears,,5]))
}

obsDout_s2<-obsDout[is.element(obsDout$siteID,intersect(obsDout$siteID,init_set2Key$Index)),]
obsDout_s2<-obsDout_s2 %>% mutate(siteID = coalesce(
  init_set2Key$siteX[match(siteID, init_set2Key$Index)]
))
obsBout_s2<-obsBout[is.element(obsBout$siteID,intersect(obsBout$siteID,init_set2Key$Index)),]
obsBout_s2<-obsBout_s2 %>% mutate(siteID = coalesce(
  init_set2Key$siteX[match(siteID, init_set2Key$Index)]
))
obsHout_s2<-obsHout[is.element(obsHout$siteID,intersect(obsHout$siteID,init_set2Key$Index)),]
obsHout_s2<-obsHout_s2 %>% mutate(siteID = coalesce(
  init_set2Key$siteX[match(siteID, init_set2Key$Index)]
))
obsVout_s2<-obsVout[is.element(obsVout$siteID,intersect(obsVout$siteID,init_set2Key$Index)),]
obsVout_s2<-obsVout_s2 %>% mutate(siteID = coalesce(
  init_set2Key$siteX[match(siteID, init_set2Key$Index)]
))
obsNout_s2<-obsNout[is.element(obsNout$siteID,intersect(obsNout$siteID,init_set2Key$Index)),]
obsNout_s2<-obsNout_s2 %>% mutate(siteID = coalesce(
  init_set2Key$siteX[match(siteID, init_set2Key$Index)]
))

init_set2New<-InitMultiSite(nYearsMS = nYears,
                            siteInfo = siteInfoNew,
                            multiInitVar = multiInitVarNew,
                            multiThin = thinningNew,
                            multiNthin = nThinningNew,
                            PAR = PAR,
                            TAir = TAir,
                            VPD = VPD,
                            Precip = Precip,
                            CO2 = CO2,
                            defaultThin = 0,
                            ClCut = 0)

init_set2New<-save(init_set2New,obsDout_s2,obsBout_s2,obsHout_s2,obsVout_s2,obsNout_s2, file="inputs/init_set2New.rdata")

### Pine - init_set3
speciesFilter <- which(init_set3$multiInitVar[,1,1]==spX1)
initPi3<-as.data.frame(init_set3$multiInitVar[speciesFilter,,1])
initPi3<-initPi3 %>%
  mutate(matches = pmap(list(D, BA),
                        function(dw, ba) filter(initdataPi,
                                                MeanD_Baweighted == dw,
                                                BasalArea == ba)$index),)
for(i in 1:length(initPi3$matches)) if (length(initPi3$matches[[i]]) == 0){
  initPi3$matches[[i]]<-NA
}

initPi3$matches<-unlist(initPi3$matches)
indN<-which(is.na(initPi3$matches)==F)
init_set3KeyPi<-data.frame(Index=initPi3$matches[indN],siteX=speciesFilter[indN])

nLayers<-max(init_set3$nLayers)
nYears<-init_set3$nYears[speciesFilter[indN]]
nSites<-length(speciesFilter[indN])

multiInitVarNew <- array(NA, dim=c(length(speciesFilter[indN]), 7,nLayers),dimnames =
                           dimnames(init_set3$multiInitVar))
siteInfoNew <- matrix(0, length(speciesFilter[indN]), 12, byrow = T,
                      dimnames = dimnames(init_set3$siteInfo))
thinningNew <- array(NA,dim=c(length(speciesFilter[indN]),4,10))

nThinningNew<-c(length(speciesFilter[indN]))

for(i in 1:length(speciesFilter[indN])){
  init_set3$multiInitVar[speciesFilter[indN][i],4,1]<-dQ[initPi3$matches[indN][i]]
  multiInitVarNew[i,,]<-init_set3$multiInitVar[speciesFilter[indN][i],,]
  siteInfoNew[i,]<-init_set3$siteInfo[speciesFilter[indN][i],]
  thinningNew[i,,]<-init_set3$thinning[speciesFilter[indN][i],,]
  nThinningNew[i]<-init_set3$nThinning[speciesFilter[indN][i]]
}
weatherNew<-init_set3$weather[unique(siteInfoNew[,2]),,,]

climIDx <- sort(unique(siteInfoNew[,2]))
siteInfoNew[,2] <- match(siteInfoNew[,2],climIDx)
maxYears <- max(init_set3$nYears[speciesFilter[indN]])
CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)

for(i in 1:length(climIDx)){
  ij <- climIDx[i]
  PAR[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,1]))
  TAir[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,2]))
  Precip[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,4]))
  VPD[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,3]))
  CO2[i,] <- as.vector(t(init_set3$weather[ij,1:maxYears,,5]))
}

obsDout_s3Pi<-obsDout[is.element(obsDout$siteID,intersect(obsDout$siteID,init_set3KeyPi$Index)),]
obsDout_s3Pi<-obsDout_s3Pi %>% mutate(siteID = coalesce(
  init_set3KeyPi$siteX[match(siteID, init_set3KeyPi$Index)]
))
obsBout_s3Pi<-obsBout[is.element(obsBout$siteID,intersect(obsBout$siteID,init_set3KeyPi$Index)),]
obsBout_s3Pi<-obsBout_s3Pi %>% mutate(siteID = coalesce(
  init_set3KeyPi$siteX[match(siteID, init_set3KeyPi$Index)]
))
obsHout_s3Pi<-obsHout[is.element(obsHout$siteID,intersect(obsHout$siteID,init_set3KeyPi$Index)),]
obsHout_s3Pi<-obsHout_s3Pi %>% mutate(siteID = coalesce(
  init_set3KeyPi$siteX[match(siteID, init_set3KeyPi$Index)]
))
obsVout_s3Pi<-obsVout[is.element(obsVout$siteID,intersect(obsVout$siteID,init_set3KeyPi$Index)),]
obsVout_s3Pi<-obsVout_s3Pi %>% mutate(siteID = coalesce(
  init_set3KeyPi$siteX[match(siteID, init_set3KeyPi$Index)]
))
obsNout_s3Pi<-obsNout[is.element(obsNout$siteID,intersect(obsNout$siteID,init_set3KeyPi$Index)),]
obsNout_s3Pi<-obsNout_s3Pi %>% mutate(siteID = coalesce(
  init_set3KeyPi$siteX[match(siteID, init_set3KeyPi$Index)]
))

init_set3NewPine<-InitMultiSite(nYearsMS = nYears,
                            siteInfo = siteInfoNew,
                            multiInitVar = multiInitVarNew,
                            multiThin = thinningNew,
                            multiNthin = nThinningNew,
                            PAR = PAR,
                            TAir = TAir,
                            VPD = VPD,
                            Precip = Precip,
                            CO2 = CO2,
                            defaultThin = 0,
                            ClCut = 0)

init_set3NewPine<-save(init_set3NewPine,obsDout_s3Pi,obsBout_s3Pi,obsHout_s3Pi,obsVout_s3Pi,obsNout_s3Pi, file="inputs/init_set3NewPine.rdata")


### init_set4
