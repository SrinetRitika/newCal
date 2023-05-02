calSet <- 10 ####calibration set (group of calibration); 1 set of calibration is composed by 20 independent calibrations
iters <- 1.5e4 ###number of iterations for each independent calibration 
fX <- 2.38
calAlg <- "DEzs"
nworkers <- 5
newV=FALSE

readName <-  paste0("chains/calOut_",calSet-1,".")
saveName <-  paste0("chains/calOut_",calSet,".")

vLocal <- FALSE