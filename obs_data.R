# Read in and process observed data
#-----
obs   <- fread('input/obs.txt',header=T)



obs <- obs[order(obs[,4],obs[,2],obs[,1],obs[,3]),]

# obs <- obs[which(is.element(obs[,1],1:nSites)),]

obs_all <- obs
obs  <- obs[-which(obs[,4]==17),]

###remove LC from calibration
# obs <- obs[-which(obs[,4]==16),]

#consider start age to find position in output array
obs_age <- obs
obs_age_all <- obs_all
obs_age$Age <- as.double(obs_age$Age)


for (i in 1:nSites){
  obs_age[siteID==i,Age:=Age-inputdata[siteID==i,age]]
  obs_age_all[siteID==i,Age:=Age-inputdata[siteID==i,age]]
}

##find run time for each site
sites <- split(obs_age,obs_age$siteID)
maxYearSite <- as.vector(sapply(lapply(sites,"[",3),max))
maxYearSite <- maxYearSite + 1

nYearsMS <- obs_age[,(max(Age)+1),by=siteID]
nYearsMS <- nYearsMS[order(nYearsMS$siteID)]
nYearsMS <- nYearsMS$V1
nYearsMS[toRemove] <- 2
# # ###remove strange data
# toRemove <- c( 75,82,
#               139:147,150:153,155:156,
#               164,352:360,
#               481:491,575)
# removeData <- which(obs_age[,1] %in% toRemove)
# obs_age <- obs_age[-removeData,]

####proc sites to remove###
obs_age <- obs_age[!siteID %in% toRemove]
# obs$siteID <- mapvalues(obs$siteID,inputdata$siteID,inputdata$siteIDnew)




outdata <- cbind(obs_age[,1],obs_age[,3],obs_age[,4],obs_age[,2],obs_age[,7])
outdata[,SpeciesID:=1]
outdata_H <- outdata[nvar==11]
outdata_D <- outdata[nvar==12]
outdata_B <- outdata[nvar==13]
outdata_V <- outdata[nvar==30]

thinSites <- which(nThinning>0)
for(i in thinSites){
  dataX <- which(outdata_H[siteID==i]$Age %in% thinning[SiteN==i]$ageThinning)
  outdata_H[siteID==i][dataX,5] <- 2
  dataX <- which(outdata_D[siteID==i]$Age %in% thinning[SiteN==i]$ageThinning)
  outdata_D[siteID==i][dataX,5] <- 2
}

outdata_H <- as.matrix(outdata_H)
outdata_D <- as.matrix(outdata_D)
outdata_B <- as.matrix(outdata_B)
outdata_V <- as.matrix(outdata_V)

obs_H <- obs_age[nvar==11]$Value
obs_D <- obs_age[nvar==12]$Value
obs_B <- obs_age[nvar==13]$Value
obs_V <- obs_age[nvar==30]$Value
# obs_data <- obs[,5]







