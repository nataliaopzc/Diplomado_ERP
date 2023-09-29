rm(list=ls())

K=2400
r=0.65
p=2.5
q=1e-3
n=20
sigma=0.2

Brms=K*(1/(1+p)^(1/p))
Rms=r/p*Brms*(1-(Brms/K)^p)
Frms=Rms/Brms
F=1.5*Frms

Biom=rep(0,1,n)
Y=rep(0,1,n)
Biom[1]=K
Y[1]=F*Biom[1]*exp(rnorm(1,0,0.1))

for (t in 2:n)
{
  Biom[t]=Biom[t-1] + r/p*Biom[t-1]*(1-(Biom[t-1]/K)^p) - Y[t-1]
  Y[t]=F*Biom[t]*exp(rnorm(1,0,0.1))
}

CPUEteo=q*Biom
CPUEdat=CPUEteo*exp(rnorm(10,0,sigma))
cv=sd(log(CPUEdat)-log(CPUEteo))


#Graficos
par(mfrow = c(1, 2))
t=seq(1,n)
plot(t,Biom,type="l",lwd=2,xlab="Tiempo", ylab="Biomasa, Capturas",col="red",ylim=c(0,K),
     main=paste("F=",round(F,2)))
lines(t,Y,type="b",col="black",lwd=2)
abline(h=Brms,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
legend(10,2200, legend=c("Biomasa","Capturas"),col=c("green","black"),lty=1,bty="n")

plot(t,CPUEteo,type="l",lwd=2,xlab="Tiempo", ylab="CPUE",col="red",ylim=c(0,q*K),
     main=paste("cv=",round(cv,3)))
lines(t,CPUEdat,type="b",col="blue",lwd=2)
abline(h=Brms,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)




