# modelos logisticos------------
rm(list=ls())

r=0.4
K=1000
p=2

B1=seq(1,50)
B2=seq(1,50)
B1[1]=10
B2[1]=10

n=50

for (t in 2:n)
{
  B1[t]=B1[t-1]+r*B1[t-1]*(1-B1[t-1]/K)
  
  B2[t]=B2[t-1]+r/p*B2[t-1]*(1-(B2[t-1]/K)^p)
}

x=seq(1,n)
plot(x,B1,type="l",lwd=2,xlab="Tiempo", ylab="Biomasa",col="red")
lines(x,B2,type="l",lwd=2,col="green")
legend("topleft",c("P-T p=1 (Schaefer)","P-T p=2.0"),col=c("red","green"),
       lwd=2,cex=1.0,bty = "n")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
