rm(list=ls()) # erasure all objects

source('read.admb.R')
source('read.admbFit.R')

 Bio=matrix(NA,6,26)
 Rec=matrix(NA,6,26)
 Fmort=matrix(NA,6,26)
 deple=matrix(NA,6,26)

 Año=seq(1993,2018)
 
for (i in 1:6)
 {
  
  
  system(paste('mestockt -ind S',i,'.dat -nox -nohess',sep=""))
  shell(paste("copy For_R.rep S",i,".rep",sep=""))
  shell(paste("copy mestock.par S",i,".par",sep=""))
  shell(paste("copy mmestock.std S",i,".std",sep=""))
  shell(paste("copy mestock.cor S",i,".cor",sep=""))
  
  
  data <-read.rep(paste("S",i,".rep",sep=""))
  attach(data)

  Bio[i,1:length(Biomasa_desovante)]=Biomasa_desovante[1,]
  Rec[i,1:length(Biomasa_desovante)]=Reclutamientos[2,]
  Fmort[i,1:length(Biomasa_desovante)]=Mort_F[1,]
  deple[i,1:length(iomasa_desovante)]=SPR[1,]
  
}
 

 Tfin=length(Biomasa_desovante[1,])
 bias=rep(NA,5)
 par(mfrow = c(2, 2))
 
 for (i in 1:5)
 {
 bias[i]=mean((Bio[i+1,1:(Tfin-(i+1))]-Bio[1,1:(Tfin-(i+1))])/Bio[1,1:(Tfin-(i+1))])
  }
 
 rho=mean(bias)
 
 matplot(Año,t(Bio)/1000,type="l",ylab="Biomasa (miles t)",lty=1,ylim=c(0,1.1/1000*max(Bio[1,])),
         main=paste("rho=",round(rho,3)))
 
 for (i in 1:5)
 {
   bias[i]=mean((Rec[i+1,1:(Tfin-(i+1))]-Rec[1,1:(Tfin-(i+1))])/Rec[1,1:(Tfin-(i+1))])
 }
 rho=mean(bias)
 matplot(Año,t(Rec),type="l",ylab="Reclutas",lty=1,ylim=c(0,max(Rec[1,])*1.1),
         main=paste("rho=",round(rho,3)))
 
 
 for (i in 1:5)
 {
   bias[i]=mean((Fmort[i+1,1:(Tfin-(i+1))]-Fmort[1,1:(Tfin-(i+1))])/Fmort[1,1:(Tfin-(i+1))])
 }
 rho=mean(bias)
 matplot(Año,t(Fmort),type="l",ylab="F ",lty=1,ylim=c(0,max(Fmort[1,])*1.1),
         main=paste("rho=",round(rho,3)))

 
 
 for (i in 1:5)
 {
   bias[i]=mean((deple[i+1,1:(Tfin-(i+1))]-deple[1,1:(Tfin-(i+1))])/deple[1,1:(Tfin-(i+1))])
 }
 rho=mean(bias)
 matplot(Año,t(deple),type="l", ylab="B/B0",lty=1,ylim=c(0,1),
         main=paste("rho=",round(rho,3)))
 abline(h=0.4,lty=2)
 
 