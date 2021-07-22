####settings


###split data for parallelization
set1 <- 1:300
set2 <- 301:600
set3 <- 601:922
set4 <- 923:948

###process parameters
parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)
param_all <- read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/parameters.csv'))
parmod <- param_all[,2] 
parmin <- param_all[,3] 
parmax <- param_all[,4] 
parnam <- as.character(param_all[,1])
npar <- length(parmod)
###error parameters
a_Hind <- which(parnam=="a_H")
b_Hind <- which(parnam=="b_H")
a_Bind <- which(parnam=="a_B")
b_Bind <- which(parnam=="b_B")
a_Dind <- which(parnam=="a_D")
b_Dind <- which(parnam=="b_D")
a_Vind <- which(parnam=="a_V")
b_Vind <- which(parnam=="b_V")
a_Hcind <- which(parnam=="a_Hc")
b_Hcind <- which(parnam=="b_Hc")
a_Wf1ind <- which(parnam=="a_Wf1")
b_Wf1ind <- which(parnam=="b_Wf1")
a_Wf2ind <- which(parnam=="a_Wf2")
b_Wf2ind <- which(parnam=="b_Wf2")
a_Acind <- which(parnam=="a_Ac")
b_Acind <- which(parnam=="b_Ac")


nCores <- ifelse(vLocal,1,4)

sets <- 1:4
likelihoods <- list()
likelihoods[[1]] <- likelihood1
likelihoods[[2]] <- likelihood2
likelihoods[[3]] <- likelihood3
likelihoods[[4]] <- likelihood4
