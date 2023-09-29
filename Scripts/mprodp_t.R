rm(list=ls())
options(warn=-1)


#Leo el archivo de datos--------------------
inputs=read.csv("data_sim8.csv",sep=",",encoding = "latin1")
attach(inputs)

# Defino los parámetros y opciones iniciales
K0=sum(Y)*.8
r0=0.5
sigma=0.15
p=1
opt_p=0 #=1 si p se estima , =0 si p es fijado (no se estima)  



#Carga funciones y librerias--------------------
library(MASS)
source("ferror_bdin.r")
source("fdinam_bdin.r")
source("fvcov_bdin.r")

parini=log(c(K0,r0,sigma,p))
datos=data.frame(Año,Y,CPUEdat,p,opt_p) # datos 
pars_fin=optim(par=parini,fn=ferror_bdin, data=datos, method="BFGS")
attach(pars_fin)

#Variables de interes--------------------
outputs=fdinam_bdin(datos,par)
attach(outputs)
K=exp(par[1])
r=exp(par[2])
sigma=exp(par[3])
p=exp(par[4])
Bmsy=K*(1/(1+p))^(1/p)
MSY=r/p*Bmsy*(1-(Bmsy/K)^p)
Fmsy=MSY/Bmsy
LL=value


#Graficos--------------------------
id=which(CPUEdat>0)
resid=log(CPUEdat[id])-log(CPUEpred[id])
resid=resid/sd(resid)

par(mfrow = c(2, 2))

id=which(CPUEdat>0)

plot(Año[id],CPUEdat[id],main=paste("Model fit  (K=",round(K,0)," r=",round(r,3),")"),ylim = c(0,max(CPUEdat[id])*1.01),
     ylab="Abundance index",xlab="Year",type="b", cex=1)#,pch = 20)
lines(Año,CPUEpred,col="red",lwd=2)
suave = smooth.spline(Año[id], CPUEdat[id], spar=0.5)
lines(suave,col="green",lwd=2,lty=1)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

barplot(resid[id]~Año[id],main="Residuales std",col="gray",ylab="Residual",xlab="Year")
abline(h = 0, col = "black")
box()
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(CPUEpred[id],resid[id],xlab="Predicted values",ylab="Residuales std",pch = 1,cex=1)#,col="gray")
abline(h=0)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(CPUEdat[id],CPUEpred[id],type="p",main=paste("r=",round(cor(log(CPUEpred[id]),log(CPUEdat[id])),2)),
     xlab="Observed",ylab="Predicted")
lines(CPUEpred,CPUEpred,type="l",col="green",lwd=2)

#Fig2-----------------------------------------------------------------------------
par(mfrow = c(2, 2))
plot(Año,Biom,ylim = c(0,max(Biom)*1.01),type="l",ylab="Biomass",xlab="Year",
     main="Biomass",pch = 16,cex=1,lwd=2)
text(min(Año)+3,Bmsy*1.2,paste("Bmsy=",round(Bmsy,0)),col="red",cex=1)
abline(h = Bmsy, col = "red",lty = 2,lwd=1)

plot(Año,Y,ylim = c(0,max(Y)*1.01),type="l",ylab="Catch",xlab="Year",
     main="Catch",pch = 16,cex=1,lwd=2)
lines(Año,G,ylim = c(0,max(G)*1.01),type="l", col="green",pch = 16,cex=1,lwd=2)
legend("topright",c("Catch","G"),col=c("black","green"),lty=1,bty="n",lwd=2)
text(min(Año)+3,MSY*1.2,paste("MSY=",round(MSY,0)),col="red",cex=1)
abline(h = MSY, col = "red",lty = 2,lwd=1)

plot(Año,Fmort,ylim = c(0,max(Fmort)*1.01),type="l",ylab="Fishing mortality",xlab="Year",
     main="Fishing mortality",pch = 16,cex=1,lwd=2)
text(min(Año)+3,Fmsy*1.2,paste("Fmsy=",round(Fmsy,3)),col="red",cex=1)
abline(h = Fmsy, col = "red",lty = 2,lwd=1)


x=seq(0,K,K/50)
y=r/p*x*(1-(x/K)^p)
plot(x,y,type="l",lwd=2,ylab="Catch",xlab="Biomass",
     main="Production", ylim=c(0,1.2*MSY))
text(Bmsy*1.1,MSY*1.1,paste("MSY=",round(MSY,0)),col="red",cex=1)
abline(v=Bmsy,col = "red",lty = 2,lwd=1)
abline(h=MSY,col = "red",lty = 2,lwd=1)

#Kobe----------------------------------------------------------------
par(mfrow = c(1, 1))
nyrs=length(Año)
SPR=Biom/K
Mort_F=Fmort
target=Bmsy/K

plot(SPR/target,Mort_F/Fmsy,pch = 16,ylab="F/Fmsy",xlab="B/Bmsy", xlim = c(0,max(SPR/target)), ylim = c(0,max(Mort_F/Fmsy)*1.5), 
     type="l",col="black",lty="dashed",main=paste("B/Bmsy=",round(SPR[nyrs]/target,2),
                                                  " F/Fmsy=",round(Mort_F[nyrs]/Fmsy,2)))
polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
polygon(c(1,1.1*max(SPR/target),1.1*max(SPR/target),1),c(0,0,1,1),col="green") #verde
polygon(c(1,1.1*max(SPR/target),1.1*max(SPR/target),1),c(1,1,1.5*max(Mort_F/Fmsy),1.5*max(Mort_F/Fmsy)),col="yellow1") #amarillo
polygon(c(0,1,1,0),c(1,1,1.5*max(Mort_F/Fmsy),1.5*max(Mort_F/Fmsy)),col="tomato1") #rojo

lines(SPR/target,Mort_F/Fmsy,pch = 16, type="b",col="black",lty="dashed",)
lines(SPR[nyrs]/target,Mort_F[nyrs]/Fmsy,type="p",col="blue",pch = 16,cex=2)
text(SPR/target*.95,Mort_F/Fmsy,paste(Año),cex=0.8)

#Variables de salida----------------
parametros=data.frame(K,r,sigma,p,Bmsy,MSY,Fmsy,LL)

variables=data.frame(Año,CPUEobs=round(CPUEdat,2),CPUEpred=round(CPUEpred,2),Capturas=Y,Biomasa=round(Biom,0),
                     F=round(Fmort,2),Produccion=round(G,0),B_Bmsy=round(Biom/Bmsy,2),F_Fmsy=round(Fmort/Fmsy,2))


vcorrel=fvcov_bdin(datos,par)
cv_par=sqrt(diag(vcorrel$vcov))
sd_par=cv_par*exp(par)

SD_par=matrix(NA,2,length(par))
SD_par[1,]=exp(par)
SD_par[2,]=sd_par
colnames(SD_par)=c('K','r','sigma','p')
rownames(SD_par)=c('est','sd(est)')

write.csv(variables, "Variables_modelo.csv")
write.csv(parametros, "Parametros.csv")
write.csv(SD_par, "sd_Parametros.csv")

muestra<-mvrnorm(1000,par,vcorrel$vcov)
par_K=exp(muestra[,1])
par_r=exp(muestra[,2])
par_sigma=exp(muestra[,3])
par_p=exp(muestra[,4])

library(psych) 
pairs(par_K~par_r+par_sigma+par_p,col="blue",main="Multiple correlations")

par(mfrow = c(2, 2))
hist(par_K,20,xlab="K",main="")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()
hist(par_r,20,xlab="r",main="")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()
hist(par_sigma,20,xlab="sigma",main="")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()
hist(par_p,20,xlab="p",main="")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()


parametros
SD_par
vcorrel
variables