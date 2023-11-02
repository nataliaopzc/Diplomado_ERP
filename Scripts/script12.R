rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=runif(1,0.2,0.8) #0.6*exp(rnorm(1,0,0.1))
Tmax=round(-log(0.01)/M)
A50=0.4*Tmax 
A95=A50+1

A50c=A50-2 
A95c=A50c+1

A50m=A50-1
A95m=A50m+1


R0=1000
h=0.75 #*exp(rnorm(1,0,0.1))
dts=0.8
Loo=150
k=0.5*M
t0=-0.5
aw=0.001
bw=3.0
Fref=M
nyears=20
sigmaR=0.4
sigmaF=0.3
sigma_cpue=0.15
sigma_Bacus=0.20
Talla=seq(30,round(1.2*Loo),4)

# vectores biologicos-pesqueros
edad=seq(1,Tmax)
Ledad=Loo*(1-exp(-k*(edad-t0)))
sd_edad=Ledad*runif(1,0.08,0.15)

wedad=aw*Ledad^bw
Om=1/(1+exp(-log(19)*(edad-A50m)/(A95m-A50m)))
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50)))
Sel_acus=1/(1+exp(-log(19)*(edad-A50c-1)/(A95c-A50c)))


# Matriz clave talla-edad inversa

dtalla=Talla[2]-Talla[1]
pdf=matrix(NA,Tmax,length(Talla))
xs=0.5*dtalla

for(i in 1: Tmax) #loop over ages
{
  for(j in 1: length(Talla)) #loop over ages
  {
    z1=((Talla[j]-xs)-Ledad[i])/sd_edad[i]
    z2=((Talla[j]+xs)-Ledad[i])/sd_edad[i]
    pdf[i,j]=pnorm(z2)-pnorm(z1)}
  pdf[i,]=pdf[i,]/sum(pdf[i,])
  
}



Ncru=matrix(NA,nyears,Tmax)
N=matrix(NA,nyears,Tmax)
F=matrix(NA,nyears,Tmax)
Z=matrix(NA,nyears,Tmax)
C=matrix(NA,nyears,Tmax)
Y=rep(NA,nyears)
B=rep(NA,nyears)
CPUE=rep(0,nyears)
Reclu=rep(NA,nyears)
Fcr=rep(NA,nyears)
Emed=rep(NA,nyears)
Bacus=rep(0,nyears)
Lmed_f=rep(0,nyears)
Lmed_n=rep(0,nyears)

yrs=seq(1995,1995+nyears-1)

# dinámica para condiciones iniciales
Reclu[1]=R0*exp(rnorm(1,0,sigmaR))
N[1,1]=Reclu[1]

for (a in 2:Tmax){
  N[1,a]=N[1,a-1]*exp(-M)
}

N[1,Tmax]=N[1,Tmax]/(1-exp(-M))


# dinámica y proyeccion

# Mortalidad por pesca anual

tlim=2000
s1=2
s2=20
  
Fcr=exp(-(yrs-tlim)^2/s1^2)
id=which(yrs>=tlim)
Fcr[id]=exp(-(yrs[id]-tlim)^2/s2^2)
Fcr=Fcr*exp(rnorm(20,0,sigmaF))

# Año 1
F[1,]=Fcr[1]*Sel
Z[1,]=F[1,]+M
C[1,]=F[1,]/Z[1,]*N[1,]*(1-exp(-Z[1,]))
Y[1]=sum(C[1,]*wedad)
B[1]=sum(N[1,]*wedad*Om*exp(-Z[1,]*dts))
CPUE[1]=0.3*sum(N[1,]*wedad*Sel)
Emed[1]=sum(C[1,]*edad)/sum(C[1,])
Ncru[1,]=N[1,]*Sel_acus*exp(-Z[1,]*dts)


alfa=4*h*R0/(5*h-1)
beta=(1-h)*B[1]/(5*h-1)



for (y in 2:nyears){
  eps=rnorm(1,0,0.15)
  Sel=1/(1+exp(-log(19)*(edad-A50*exp(eps))/((A95-A50)*exp(eps))))
  F[y,]=Fcr[y]*Sel
  Z[y,]=F[y,]+M
  
  Reclu[y]=alfa*B[y-1]/(beta+B[y-1])*exp(rnorm(1,0,sigmaR))
  N[y,1]=Reclu[y]

  for (a in 2:Tmax){
  N[y,a]=N[y-1,a-1]*exp(-Z[y-1,a-1])
  }
  
  N[y,Tmax]=N[y,Tmax]+N[y-1,Tmax-1]*exp(-Z[y-1,Tmax-1])

  C[y,]=F[y,]/Z[y,]*N[y,]*(1-exp(-Z[y,]))
  Emed[y]=sum(C[y,]*edad)/sum(C[y,])
  
  Y[y]=sum(C[y,]*wedad)
  B[y]=sum(N[y,]*wedad*Om*exp(-Z[y,]*dts))
  Bacus[y]=0.9*sum(N[y,]*wedad*Sel_acus*exp(-Z[y,]*dts))*exp(rnorm(1,0,sigma_Bacus))
  CPUE[y]=0.3*sum(N[y,]*wedad*Sel)*exp(rnorm(1,0,sigma_cpue))
  Ncru[y,]=N[y,]*Sel_acus*exp(-Z[y,]*dts)

}




nmf=5000
nmn=1000

Ntalla=t(pdf)%*%t(Ncru)
Ctalla=t(pdf)%*%t(C)

propl_f=Ctalla
propl_n=Ntalla

#par(mfrow = c(5, 4))

 for (y in 1:nyears){
 
 propl_f[,y]=rmultinom(1, nmf, Ctalla[,y]/sum(Ctalla[,y]))/nmf # muestreo con error multinomial
 propl_n[,y]=rmultinom(1, nmn, Ntalla[,y]/sum(Ntalla[,y]))/nmn # muestreo con error multinomial
 Lmed_f[y]= sum(propl_f[,y]*Talla)/sum(propl_f[,y])
 Lmed_n[y]= sum(propl_n[,y]*Talla)/sum(propl_n[,y])
 
 
 #plot(Talla,propl_f[,y],type="l",main=paste(y))
 #lines(Talla,propl_n[,y],col="green")
 
 }

  
  

Bacus[1:10]=0
Ctalla_flo= propl_f;Ctalla_flo[,1:5]=0
Ntalla_acu= propl_n; Ntalla_acu[,1:10]=0; Lmed_n[1:10]=0

Capt_l=round(Ctalla_flo*1500)
Cru_l=round(Ntalla_acu*1500)

indices=data.frame(Yrs=yrs,Catch=Y,CPUE=CPUE/max(CPUE),Bcru=Bacus)


matCap=data.frame(t(Capt_l))
rownames(matCap)=yrs
colnames(matCap)=Talla

matNcru=data.frame(t(Cru_l))
rownames(matNcru)=yrs
colnames(matNcru)=Talla


params=c(106.1,117.5,log(aw),bw,dts,0.5,dts,Loo,k,M)
tabla1 <- matrix(ncol=1, round(params,3))
rownames(tabla1) <- c("L50m", "L95m","log_aw","bw","dt_spawn", "dt_cpue", "dt_crucero","Loo","k","M" )
colnames(tabla1)<-"Value"

tabla2=data.frame(yrs=yrs,Biom=B,Reclu=Reclu,mortF=Fcr)


write.csv(tabla1, 'Parsbiol.csv',row.names = T)

write.csv(tabla2, 'Vars_vd.csv',row.names = T)

write.csv(indices, 'indices.csv',row.names = T)

write.csv(matCap, 'MatCap.csv', row.names = T)

write.csv(matNcru, 'MatNcru.csv', row.names = T)



par(mfrow = c(2, 3))

barplot(Reclu~yrs, main="Reclutamiento")
box()
abline(h=R0,lty=2)

plot(yrs,Fcr,type="l",col="red",lwd=2,main="Mort por pesca", ylim=c(0,max(Fcr)*1.1))
abline(h=Fref,lty=2)

barplot(Y~yrs,col="lightblue",main="Capturas", ylim=c(0,max(Y)*1.1))
box()
abline(h=mean(Y),lty=2)

maxB=max(Bacus[Bacus>0])
Bacus[Bacus==0]=NaN
plot(yrs,B/max(B),type="l",col="green",lwd=2,ylim=c(0,max(B/B[1])*1.1),main="Biomasa y CPUE")
lines(yrs,CPUE/max(CPUE))
lines(yrs,Bacus/maxB,col="red",type="b")
abline(h=0.4,lty=2)

Lmed_n[Lmed_n==0]=NaN
plot(yrs,Lmed_f,type="l",lwd=2,main="Talla promedio", ylab="Talla", ylim=c(min(Lmed_n[seq(nyears-9,nyears)]),max(Lmed_f)))
lines(yrs,Lmed_n,type="l",lwd=2,col="green")

plot(Talla,propl_f[,nyears],type="l",ylab="Frecuencia",main="Frec tallas",col="black",lwd=2)
lines(Talla,propl_f[,1], col="green",lwd=2)
legend('topleft',c('final','inicial'),col=c('black','green'),lty=1,lwd=2,bty="n")



