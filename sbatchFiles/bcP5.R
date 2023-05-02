setwd("/scratch/project_2000994/calibrations/srinet/newCal/sbatchFiles/")
source("genSet.r")
pNum <- "5"
fX <- fX*3
###load data
lastCal <- paste0(readName,pNum,".rdata")
###writeName
newCal <- paste0(saveName,pNum,".rdata")

# if(vLocal){
#   nworkers <- 4 # required cores
# }else{
#   nworkers <- 30 # required cores
#   iters <- 6e2
#   #Make a cluster with nworkers
#   cl <- parallel::makeCluster(nworkers)
# }

setwd("/scratch/project_2000994/calibrations/srinet/newCal/")

startX <- Sys.time()
source("Rsrc/2.1_BC_model.R")
endX <- Sys.time()
timeX = endX- startX
print(timeX)