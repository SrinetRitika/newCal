library(ggplot2)
library(Rprebasso)

## vapu spruce.
vapu_S<-read.csv("inputs/VAPU_spruce.csv")
nData_S <- length(vapu_S$plotNo)
pars <- pCROB[c(11,15,38),2] #z,rhof2,ksi
# pars<-c(1.927,408.92,0.08178)
Wf <- matrix(0.,nData_S,2)
inputs <- matrix(0.,nData_S,3,dimnames = list(NULL,c("BA","H","Hc")))
inputs[,1] <- sample(vapu_S$BA,nData_S)
inputs[,2] <- sample(vapu_S$H.m,nData_S)
inputs[,3] <- inputs[,2]-0.3*inputs[,2]
As <- matrix(0., nData_S,2)
test <- .Fortran("calWf",pars=as.double(pars),
                 Wf=as.matrix(Wf),
                 inputs=as.matrix(inputs),
                 nData=as.integer(nData_S),
                 As=as.matrix(As))

Lc <- inputs[,2]-inputs[,3]
Wf_Lc <- .Fortran("calWf_fLc",pars=as.double(pars[c(3,1)]),
                 Wf=as.double(Wf),
                 nData=as.integer(nData_S),
                 Lc=as.double(Lc)
                 )

As <- vapu_S$Ac
Wf_As <- .Fortran("calWf_fA",pars=as.double(pars[2]),
                 Wf=as.double(Wf),
                 nData=as.integer(nData_S),
                 As=as.double(As)
)

# As <- vapu_S$Ac
As_Lc <- .Fortran("calAs_fLc",pars=as.double(pars[c(3,1,2)]),
                  As=as.double(As),
                  nData=as.integer(nData_S),
                  Lc=as.double(Lc)
)

calAs_fLc(pars,As,nData,Lc)


data<-as.data.frame(cbind(test$As[,1],test$As[,2],vapu_S$Ac,test$Wf[,1],test$Wf[,2],vapu_S$Wf.1))
colnames(data)<-c("As1","As2","As_m","Wf1","Wf2","Wf_m")
p1<-ggplot(data=data,aes(x=As_m,y=As1))+geom_point()+geom_abline()+geom_point(data=data,aes(x=As_m,y=As2),col="red")
p1
p2<-ggplot(data=data,aes(x=Wf_m,y=Wf1))+geom_point()+geom_abline()+geom_point(data=data,aes(x=Wf_m,y=Wf2),col="red")
p2
## vapu pine.
vapu_P<-read.csv("inputs/VAPU_pine.csv")
nData_P <- length(vapu_P$plot)
pars <- pCROB[c(11,15,38),1] #z,rhof2,ksi
pars<-c(1.9186,211.29,0.084516)
Wf <- matrix(0.,nData_P,2)
inputs <- matrix(0.,nData_P,3,dimnames = list(NULL,c("BA","H","Hc")))
inputs[,1] <- sample(vapu_P$BA,nData_P)
inputs[,2] <- sample(vapu_P$H,nData_P)
inputs[,3] <- inputs[,2]-0.3*inputs[,2]
As <- matrix(0., nData_P,2)
test <- .Fortran("calWf",pars=as.double(pars),
                 Wf=as.matrix(Wf),
                 inputs=as.matrix(inputs),
                 nData=as.integer(nData_P),
                 As=as.matrix(As))
data_p<-as.data.frame(cbind(test$As[,1],test$As[,2],vapu_P$Ac,test$Wf[,1],test$Wf[,2],vapu_P$Wf.1))
colnames(data_p)<-c("As1","As2","As_m","Wf1","Wf2","Wf_m")
p1<-ggplot(data=data_p,aes(x=As_m,y=As1))+geom_point()+geom_abline()+geom_point(data=data_p,aes(x=As_m,y=As2),col="red")
p1
p2<-ggplot(data=data_p,aes(x=Wf_m,y=Wf1))+geom_point()+geom_abline()+geom_point(data=data_p,aes(x=Wf_m,y=Wf2),col="red")
p2