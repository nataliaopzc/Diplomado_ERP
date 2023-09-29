rm(list=ls())

r=0.45
K=1000
p=1
n=50

B1=seq(1,n)
Capt=seq(1,n)

B1[1]=K
mu=0.225


Capt[1]=mu*B1[1]

for (t in 2:n)
{
  
  B1[t]=(B1[t-1]+r/p*B1[t-1]*(1-(B1[t-1]/K)^p)-Capt[t-1])*exp(rnorm(1,0,0.0))
  
  Capt[t]=mu*B1[t]
  

}

x=seq(1,n)
plot(x,B1,type="l",lwd=2,xlab="Tiempo", ylab="Biomasa",col="red",ylim=c(0,K),
     main=paste("mu=",mu))
lines(x,Capt,lwd=2)
abline(h=K/2,col="red",lty=2)
abline(h=K*r/4,col="black",lty=2)

legend("top",c("Biomasa","Captura"),col=c("red","black"),lwd=2,cex=1.0,bty = "n")