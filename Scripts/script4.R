# Curva de producción para distintos niveles de productividad
rm(list=ls())

r=0.45
K=1000

B=seq(0,K,10) # vector de valores ceros entre cero y K con paso de 10

# Curvas de producción
p=1e-5
G1=r/p*B*(1-(B/K)^p)
Brms1=K*(1/(p+1))^(1/p)
rms1=r/p*Brms1*(1-(Brms1/K)^p)

p=1
G2=r/p*B*(1-(B/K)^p)
Brms2=K*(1/(p+1))^(1/p)
rms2=r/p*Brms2*(1-(Brms2/K)^p)

p=3
G3=r/p*B*(1-(B/K)^p)
Brms3=K*(1/(p+1))^(1/p)
rms3=r/p*Brms3*(1-(Brms3/K)^p)


# Graficos
plot(B,G2,col="green",lwd=2,type="l",
     ylab="Producción",xlab="Biomasa", main=paste("K= ",K,"r= ",r))
abline(v=Brms1,lty=2,col="green")
abline(h=rms1,lty=2,col="green")

lines(B,G2,col="red",lwd=2)
abline(v=Brms2,lty=2,col="red")
abline(h=rms2,lty=2,col="red")

lines(B,G3,col="black",lwd=2)-
abline(v=Brms3,lty=2)
abline(h=rms3,lty=2)

legend(x=800,y=150,c("p=1e-5","p=1","p=3"),
       col=c("green","red","black"),lwd=2,bty="n")

