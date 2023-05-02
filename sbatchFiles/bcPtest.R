setwd("/scratch/project_2000994/calibrations/srinet/newCal/sbatchFiles/")
source("genSet.r")
pNum <- "1"
###load data
lastCal <- paste0(readName,pNum,".rdata")
###writeName
newCal <- paste0(saveName,pNum,".rdata")

# if(vLocal){
#   nworkers <- 4 # required cores
# }else{
#   nworkers <- 30 # required cores
#   #Make a cluster with nworkers
#   cl <- parallel::makeCluster(nworkers)
# }

setwd("/scratch/project_2000994/calibrations/srinet/newCal/")

startX <- Sys.time()
source("Rsrc/test_BC_model.R")
endX <- Sys.time()
timeX = endX- startX
print(timeX)