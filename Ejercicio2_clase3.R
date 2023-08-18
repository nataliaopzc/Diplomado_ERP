## Escenario donde agregamos error de observacón
rm(list=ls())

# Defino los parámetros poblacionales
K=2400
r=0.65
p=3
q=1e-3
sigma=0.15 # error de observación
sigmay=0.2
n=20 #número de años

# Defino los PBR
Brms=K*(1/(p+1))^(1/p)
RMS=(r/p)*Brms*(1-(Brms/K)^p)
Frms=RMS/Brms

# Defino las condiciones iniciales 
Biom=rep(0,1,n) #repetir los valores 0 de 1 hasta n
Y=rep(0,1,n)
Biom[1]=K
F=1.5*Frms
Y[1]=0 #Pesquería en condiciones virginales #F*Biom[1]*exp(rnorm(n,0,sigmay))

# Proyecto la biomasa y las capturas desde el año 2
for (t in 2:n)
{
  Biom[t]=Biom[t-1] + (r/p)*Biom[t-1]*(1-(Biom[t-1]/K)^p) - Y[t-1]
  Y[t]=1.5*Frms*Biom[t]*exp(rnorm(n,0,sigmay))
}

CPUEteo=q*Biom
CPUEdat=CPUEteo*exp(rnorm(n,0,sigma))
CPUErms=q*Brms

cv=sd(log(CPUEdat)-log(CPUEteo))
#print(paste('cv=',round(cv,3)))

# Gráficos
par(mfrow= c(1,2))

t=seq(1,n)
plot(t,Biom,type="l", col='red',lwd='2',ylim=c(0,K))
lines(t,Y,type="b", col='black',lwd='2')
abline(h=Brms,col='green',lwd='2',lty=2)

plot(t,CPUEteo,type="l", col='red',lwd='2',ylim=c(0,max(CPUEteo)),
     main=paste("cv=",round(cv,3)))
lines(t,CPUEdat,type="b", col='black',lwd='2')
abline(h=CPUErms,col='green',lwd='2',lty=2)
