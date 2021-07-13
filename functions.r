##I have intergrated the old growth data into initprebas and made it data set4
##questions:how to set pvalues in the likelihood function? a serial of numbers for each parameter

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
###init_set4
init_set4$pCROBAS <- pCROB
init_set4$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
init_set4$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
init_set4$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]

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
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB*3+1)]+pValues[(nparCROB*3+2)]*output[Hdata_s1$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB*3+3)]+pValues[(nparCROB*3+4)]*output[Ddata_s1$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB*3+5)]+pValues[(nparCROB*3+6)]*output[Bdata_s1$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB*3+7)]+pValues[(nparCROB*3+8)]*output[Hcdata_s1$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB*3+9)]+pValues[(nparCROB*3+10)]*output[Vdata_s1$outData]))
  
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
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB*3+1)]+pValues[(nparCROB*3+2)]*output[Hdata_s2$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB*3+3)]+pValues[(nparCROB*3+4)]*output[Ddata_s2$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB*3+5)]+pValues[(nparCROB*3+6)]*output[Bdata_s2$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB*3+7)]+pValues[(nparCROB*3+8)]*output[Hcdata_s2$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB*3+9)]+pValues[(nparCROB*3+10)]*output[Vdata_s2$outData]))
  
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
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB*3+1)]+pValues[(nparCROB*3+2)]*output[Hdata_s3$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB*3+3)]+pValues[(nparCROB*3+4)]*output[Ddata_s3$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB*3+5)]+pValues[(nparCROB*3+6)]*output[Bdata_s3$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB*3+7)]+pValues[(nparCROB*3+8)]*output[Hcdata_s3$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB*3+9)]+pValues[(nparCROB*3+10)]*output[Vdata_s3$outData]))
  
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

##########old growth likelihood
likelihood4 <- function(pValues){
  ###init_set4
  init_set4$pCROBAS <- pCROB
  init_set4$pCROBAS[parSel,1] <- pMAP[1:nparCROB]
  init_set4$pCROBAS[parSel,2] <- pMAP[(nparCROB + 1):(nparCROB*2)]
  init_set4$pCROBAS[parSel,3] <- pMAP[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set4)$multiOut
  
  outdata_B2 <- Bdata_s4$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s4$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s4$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s4$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s4$outData]-Hdata_s4$obs
  diff_Hc <- output[Hcdata_s4$outData]-Hcdata_s4$obs
  diff_D <- output[Ddata_s4$outData]-Ddata_s4$obs
  diff_B <- out_B-Bdata_s4$obs
  diff_V <- out_V-Vdata_s4$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[(nparCROB*3+1)]+pValues[(nparCROB*3+2)]*output[Hdata_s4$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[(nparCROB*3+3)]+pValues[(nparCROB*3+4)]*output[Ddata_s4$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[(nparCROB*3+5)]+pValues[(nparCROB*3+6)]*output[Bdata_s4$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[(nparCROB*3+7)]+pValues[(nparCROB*3+8)]*output[Hcdata_s4$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[(nparCROB*3+9)]+pValues[(nparCROB*3+10)]*output[Vdata_s4$outData]))
  
  ##vapu data
  As_p_obs <- vapu_P$Ac
  wf_p_obs <- vapu_P$Wf.1 ## wf.1 is the carbon, wf is the kg.
  Lc_p <- vapu_P$Hc
  As_s_obs <- vapu_S$Ac
  wf_s_obs <- vapu_S$Wf.1
  Lc_s<-vapu_S$Hc.m
  
  ksi_p <- init_set4$pCROBAS[38,1]
  rhof_p <- init_set4$pCROBAS[15,1]
  z_p <- init_set4$pCROBAS[11,1]
  Wf1_p <- rhof_p*As_p_obs
  Wf2_p <- ksi_p*Lc_p^z_p
  As_p <- ksi_p/rhof_p*Lc_p^z_p
  
  diff_wf1_p <- Wf1_p-wf_p_obs
  diff_wf2_p <- Wf2_p-wf_p_obs
  diff_As_p <- As_p-As_p_obs
  
  ll_wf1_p <- sum(Sivia_log(diff_wf1_p,sd = abs(diff_wf1_p*0.3)))
  ll_wf2_p <- sum(Sivia_log(diff_wf2_p,sd = abs(diff_wf2_p*0.3)))
  ll_As_p <- sum(Sivia_log(diff_As_p,sd = abs(diff_As_p*0.3)))
  
  ksi_s <- init_set4$pCROBAS[38,2]
  rhof_s <- init_set4$pCROBAS[15,2]
  z_s <- init_set4$pCROBAS[11,2]
  Wf1_s <- rhof_s*As_s_obs
  Wf2_s <- ksi_s*Lc_s^z_s
  As_s <- ksi_s/rhof_s*Lc_s^z_s
  
  diff_wf1_s <- Wf1_s-wf_s_obs
  diff_wf2_s <- Wf2_s-wf_s_obs
  diff_As_s <- As_s-As_s_obs
  
  ll_wf1_s <- sum(Sivia_log(diff_wf1_s,sd = abs(diff_wf1_s*0.3)))
  ll_wf2_s <- sum(Sivia_log(diff_wf2_s,sd = abs(diff_wf2_s*0.3)))
  ll_As_s <- sum(Sivia_log(diff_As_s,sd = abs(diff_As_s*0.3)))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V,ll_wf1_p,ll_wf2_p,ll_As_p,ll_wf1_s,ll_wf2_s,ll_As_s)
  
  return(loglikelihood)
}

