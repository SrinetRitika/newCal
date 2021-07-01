parSel <- c(1:18,31:34,38,41)
nparCROB <- length(parSel)

initPrebas$pCROBAS <- pCROB
initPrebas$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
initPrebas$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
initPrebas$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
###init_set2
init_set2$pCROBAS <- pCROB
init_set2$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set2$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set2$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
###init_set3
init_set3$pCROBAS <- pCROB
init_set3$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set3$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set3$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

likelihood1 <- function(pValues){
  ###init_set1
  init_set1$pCROBAS <- pCROB
  init_set1$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
  init_set1$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
  init_set1$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set1)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s1$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s1$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s1$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s1$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s1$outData]-Hdata_s1$obs
  diff_Hc <- output[Hcdata_s1$outData]-Hcdata_s1$obs
  diff_D <- output[Ddata_s1$outData]-Ddata_s1$obs
  diff_B <- out_B-Bdata_s1$obs
  diff_V <- out_V-Vdata_s1$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB+1)]+pValues[(nparCROB+2)]*output[Hdata_s1$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB+3)]+pValues[(nparCROB+4)]*output[Ddata_s1$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB+5)]+pValues[(nparCROB+6)]*output[Bdata_s1$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB+7)]+pValues[(nparCROB+8)]*output[Hcdata_s1$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB+9)]+pValues[(nparCROB+10)]*output[Vdata_s1$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  
  return(loglikelihood)
}


likelihood2 <- function(pValues){
  ###init_set2
  init_set2$pCROBAS <- pCROB
  init_set2$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
  init_set2$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
  init_set2$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set2)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s2$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s2$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s2$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s2$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s2$outData]-Hdata_s2$obs
  diff_Hc <- output[Hcdata_s2$outData]-Hcdata_s2$obs
  diff_D <- output[Ddata_s2$outData]-Ddata_s2$obs
  diff_B <- out_B-Bdata_s2$obs
  diff_V <- out_V-Vdata_s2$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB+1)]+pValues[(nparCROB+2)]*output[Hdata_s2$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB+3)]+pValues[(nparCROB+4)]*output[Ddata_s2$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB+5)]+pValues[(nparCROB+6)]*output[Bdata_s2$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB+7)]+pValues[(nparCROB+8)]*output[Hcdata_s2$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB+9)]+pValues[(nparCROB+10)]*output[Vdata_s2$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  
  return(loglikelihood)
}


likelihood3 <- function(pValues){
  ###init_set3
  init_set3$pCROBAS <- pCROB
  init_set3$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
  init_set3$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
  init_set3$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set3)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s3$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s3$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s3$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s3$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s3$outData]-Hdata_s3$obs
  diff_Hc <- output[Hcdata_s3$outData]-Hcdata_s3$obs
  diff_D <- output[Ddata_s3$outData]-Ddata_s3$obs
  diff_B <- out_B-Bdata_s3$obs
  diff_V <- out_V-Vdata_s3$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB+1)]+pValues[(nparCROB+2)]*output[Hdata_s3$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB+3)]+pValues[(nparCROB+4)]*output[Ddata_s3$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB+5)]+pValues[(nparCROB+6)]*output[Bdata_s3$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB+7)]+pValues[(nparCROB+8)]*output[Hcdata_s3$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB+9)]+pValues[(nparCROB+10)]*output[Vdata_s3$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  
  return(loglikelihood)
}

