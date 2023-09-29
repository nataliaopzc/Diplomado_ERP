rm(list=ls())

r=0.65
K=2400
n=50

Biom=rep(0,1,n)
G=rep(0,1,n)

Biom[1]=K
G[1]=r*Biom[1]*(1-Biom[1]/K)

Captura=350

for (t in 2:n)
{
  Biom[t]=Biom[t-1] + G[t-1] - Captura
  G[t]=r*Biom[t]*(1-Biom[t]/K)
  
}

par(mfrow = c(1, 2))

plot(Biom,type="l",lwd=2,xlab="Tiempo", ylab="Biomasa",col="red",ylim=c(0,K),
     main=paste("Captura=",Captura," Bf/K=",round(Biom[n]/K,2)))
abline(h=0.5*K,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

B=seq(0,K,50)
Grow=r*B*(1-B/K)
plot(B,Grow,type="l",xlab="Biomasa",ylab="Producción",col="gray")
lines(Biom,G,type="b",col="red")

lines(Biom,Captura*rep(1,1,n),type="p")
abline(h=r*K/4,lwd=1,col="gray")
abline(v=K/2,lwd=1,col="gray")


