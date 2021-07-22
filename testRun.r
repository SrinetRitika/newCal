

####Run model in parallel
library(ggplot2)
# library(Metrics)
library(data.table)
library(parallel)
library(BayesianTools)
library(ggpubr)

####set project library if running on CSC
vLocal <- TRUE #### flag for runs on CSC(FALSE) or on laptop(TRUE)
if(!vLocal){
  setwd("/scratch/project_2000994/calibrations/all")
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

# setwd("/scratch/project_2000994/calibrations/newCal")
setwd("C:/Users/checcomi/Documents/github/newCal")

source("functions.r") ###run using url
source("settings.r") ###run using url

# load("inputs/initPrebas_old growth.rdata")
load("inputs/init_set1.rdata")
load("inputs/init_set2.rdata")
load("inputs/init_set3.rdata")
load("inputs/init_set4.rdata")

load("outCal/pMAP.rdata")


# modOut <- list()
# 
###Run model 
startX <- Sys.time()
likX <- likelihood(pMAP)
endX <- Sys.time()
timeX = endX- startX
print("parallel runs")
print(likX)
print(timeX)

startX <- Sys.time()
ll1 <- likelihood1(pMAP)
ll2 <- likelihood2(pMAP)
ll3 <- likelihood3(pMAP)
ll4 <- likelihood4(pMAP)
llS <- ll1 + ll2 +ll3+ll4
endX <- Sys.time()
timeX = endX- startX
print("sequential runs")
print(llS)
print(timeX)
