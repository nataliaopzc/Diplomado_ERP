rm(list=ls())

B=seq(0,1000,10)
r=0.3
K=1000
p=5
b=0.01


G1=(r+1)*B
G2=r*B*(1-B/K)
G3=r/p*B*(1-(B/K)^p)
G4=(r+1)*B/(1+b*B)

plot(B,G1/max(G1),type="l",lwd=2,ylab="Relative G", xlab="Biomass")
lines(B,G2/max(G2),col="red",lwd=2)
lines(B,G3/max(G3),col="green",lwd=2)
lines(B,G4/max(G4),col="blue",lwd=2)
legend("topleft",c("lineal", "P-T p=1 (Schaefer)","P-T p=2.0","B-H"),lty =c(1,1,1),col=c("black","red","green","blue"),
       lwd=2,cex=1.0,bty = "n")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)




