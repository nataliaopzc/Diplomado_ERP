# Curva de producción para distintos niveles de productividad
rm(list=ls())
options(warn=-1)
setwd('C:/Users/n_aty/OneDrive/Documentos/Diplomado/Clases/1° Control/')

#Leo el archivo de datos--------------------
#inputs=read.csv("dataNO.csv",sep=";",encoding = "latin1")
#attach(inputs)

r=0.3228347
K1=2732.97
#K=1000

B1=seq(0,K1,10) # vector de valores ceros entre cero y K con paso de 10

# Curvas de producción
p=1e-3
G1=r/p*B1*(1-(B1/K1)^p)
Brms1=K1*(1/(p+1))^(1/p)
rms1=r/p*Brms1*(1-(Brms1/K1)^p)

Brms1_1=0.4*K1
rms1_1=r/p*Brms1_1*(1-(Brms1_1/K1)^p)

# Graficos
plot(B1,G1,ylim = c(0, 450),col="green",lwd=2,type="l",
     ylab="Producción",xlab="Biomasa")
text(Brms1_1+3,rms1_1*0.95,paste("40%B0=",round(rms1_1,0)),col="green",cex=0.8)
abline(v=Brms1_1,lty=2,col="green")
abline(h=rms1_1,lty=2,col="green")

########################################
K2=2594.131
r=0.5661019

B2=seq(0,K2,10) 
p=1
G2=r/p*B2*(1-(B2/K2)^p)
Brms2=K2*(1/(p+1))^(1/p)
rms2=r/p*Brms2*(1-(Brms2/K2)^p)

Brms2_1=0.4*K2
rms2_1=r/p*Brms2_1*(1-(Brms2_1/K2)^p)

lines(B2,G2,col="red",lwd=2)
text(Brms2_1-200,rms2_1*1.03,paste("40%B0=",round(rms2_1,0)),col="red",cex=0.8)
abline(v=Brms2_1,lty=2,col="red")
abline(h=rms2_1,lty=2,col="red")

#########################################
K3=2479.113
r=1.081457

B3=seq(0,K3,10) 
p=3
G3=r/p*B3*(1-(B3/K3)^p)
Brms3=K3*(1/(p+1))^(1/p)
rms3=r/p*Brms3*(1-(Brms3/K3)^p)

Brms3_1=0.4*K3
rms3_1=r/p*Brms3_1*(1-(Brms3_1/K3)^p)

lines(B3,G3,col="black",lwd=2)
text(Brms3_1-250,rms3_1*1,paste("40%B0=",round(rms3_1,0)),col="black",cex=0.8)
abline(v=Brms3_1,lty=2)
abline(h=rms3_1,lty=2)
#abline(v=Brms3/,lty=2)

legend(x=2100,y=410,c("p=1e-3","p=1","p=3"),
       col=c("green","red","black"),lwd=2,bty="n")

