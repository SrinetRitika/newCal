####settings

parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)

###split data for parallelization
set1 <- 1:300
set2 <- 301:600
set3 <- 601:922
set4 <- 923:948
