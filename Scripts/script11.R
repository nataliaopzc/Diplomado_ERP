rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=0.3
Tmax=round(-log(0.01)/M)
A50=0.4*Tmax 
A95=A50+1
R0=10000
dts=0.8
Loo=150
k=0.5*M
t0=-0.5
aw=0.001
bw=3.0
A50m=0.67*Tmax
A95m=A50m+1
Fref=0.5*M
nyears=20
sigmaR=0.6
sigmaF=0.2


# vectores biologicos-pesqueros
edad=seq(1,Tmax)
Ledad=Loo*(1-exp(-k*(edad-t0)))
wedad=aw*Ledad^bw
Om=1/(1+exp(-log(19)*(edad-A50m)/(A95m-A50m)))
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50)))


N=matrix(NA,nyears,Tmax)
F=matrix(NA,nyears,Tmax)
Z=matrix(NA,nyears,Tmax)
C=matrix(NA,nyears,Tmax)
Y=rep(NA,nyears)
B=rep(NA,nyears)
CPUE=rep(NA,nyears)
Reclu=rep(NA,nyears)
Fcr=rep(NA,nyears)
Emed=rep(NA,nyears)

# dinámica para condiciones iniciales
Reclu[1]=R0
N[1,1]=Reclu[1]

for (a in 2:Tmax){
  N[1,a]=N[1,a-1]*exp(-M)
}

N[1,Tmax]=N[1,Tmax]/(1-exp(-M))


# dinámica y proyeccion

# Año 1
Fcr[1]=Fref*exp(rnorm(1,0,sigmaF))
F[1,]=Fcr[1]*Sel
Z[1,]=F[1,]+M
C[1,]=F[1,]/Z[1,]*N[1,]*(1-exp(-Z[1,]))
Y[1]=sum(C[1,]*wedad)
B[1]=sum(N[1,]*wedad*Om*exp(-Z[1,]*dts))
CPUE[1]=0.3*sum(N[1,]*wedad*Sel)
Emed[1]=sum(C[1,]*edad)/sum(C[1,])

for (y in 2:nyears){
  Fcr[y]=Fref*exp(rnorm(1,0,sigmaF))
  F[y,]=Fcr[y]*Sel
  Z[y,]=F[y,]+M
  
  Reclu[y]=R0*exp(rnorm(1,0,sigmaR))
  N[y,1]=Reclu[y]

  for (a in 2:Tmax){
  N[y,a]=N[y-1,a-1]*exp(-Z[y-1,a-1])
  }
  
  N[y,Tmax]=N[y,Tmax]+N[y-1,Tmax-1]*exp(-Z[y-1,Tmax-1])

  C[y,]=F[y,]/Z[y,]*N[y,]*(1-exp(-Z[y,]))
  Emed[y]=sum(C[y,]*edad)/sum(C[y,])
  
  Y[y]=sum(C[y,]*wedad)
  B[y]=sum(N[y,]*wedad*Om*exp(-Z[y,]*dts))
  CPUE[y]=0.3*sum(N[y,]*wedad*Sel)
  
}


yrs=seq(1,nyears)
par(mfrow = c(2, 3))

barplot(Reclu~yrs, main="Reclutamiento")
box()
abline(h=R0,lty=2)

plot(yrs,Fcr,type="l",col="red",lwd=2,main="Mort por pesca", ylim=c(0,max(Fcr)*1.1))
abline(h=Fref,lty=2)

barplot(Y~yrs,col="lightblue",main="Capturas", ylim=c(0,max(Y)*1.1))
box()
abline(h=mean(Y),lty=2)


plot(yrs,B/max(B),type="l",col="green",lwd=2,ylim=c(0,max(B/B[1])*1.1),main="Biomasa y CPUE")
lines(yrs,CPUE/max(CPUE))
abline(h=0.4,lty=2)

plot(yrs,Emed,type="l",lwd=2,main="Edad promedio", ylab="Edad")

barplot(C[2,]/sum(C[2,])~edad,ylab="Frecuencia",main="Frec edades")
barplot(C[nyears,]/sum(C[nyears,])~edad, add = T,col = "lightblue", ylim=c(0,0.3))
box()



