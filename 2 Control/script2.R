rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=0.33
Tmax=12
A50=2.4 # edad completamente reclutada # 50%
A95=3

######################
#A50_2=3.15 # edad completamente reclutada # 50%
#A95_2=4

##########################33
#Fref=0.12 #para calculat el ratio de agotamaiento (Biomasa actual /bomasa virginial (F=0)). Ahí te da el % de biomasa en el agua y puedes ver qué F te sirve

#Fref=2.5*M #para calculat el ratio de agotamaiento (Biomasa actual /bomasa virginial (F=0)). Ahí te da el % de biomasa en el agua y puedes ver qué F te sirve
R0=1
#dtf=0.5
dts=0.5833 #picos reproductivos
Loo=800.4
k=0.14
t0=-0.918
aw=0.001 #parámetro peso/talla
bw=3.0 #potencia peso/talla
A50m=2.8 #(año 2021)
A95m=3.5 # 
h=0.65 #(stock poco resiliente)


# vectores biologicos
edad=seq(1,Tmax)
Ledad=Loo*(1-exp(-k*(edad-t0)))
wedad=aw*Ledad^bw
N=rep(0,Tmax)
Om=1/(1+exp(-log(19)*(edad-A50m)/(A95m-A50m))) #modelo logístico a la edad
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50))) # selectividad

Fref=seq(0,5*M,M/10) # genera un vector de valores de F
YPR=rep(0,length(Fref))
BPR=rep(0,length(Fref))
YPR_eq=rep(0,length(Fref))
BPR_eq=rep(0,length(Fref))
Reclu=rep(0,length(Fref))


for (i in 1:length(Fref)){
  
  
  F=Fref[i]*Sel
  Z=F+M
  
  # dinámica para B0
  N[1]=R0
  
  for (a in 2:Tmax){
    N[a]=N[a-1]*exp(-M)
  }
  
  N[Tmax]=N[Tmax]/(1-exp(-M)) # número del grupo plus
  B0=sum(N*wedad*Om*exp(-M*dts))
  
  
  # dinámica en F
  
  N[1]=R0
  
  for (a in 2:Tmax){
    N[a]=N[a-1]*exp(-Z[a-1])
  }
  
  N[Tmax]=N[Tmax]/(1-exp(-Z[Tmax]))
  
  C=F/Z*N*(1-exp(-Z))
  
  YPR[i]=sum(C*wedad) #Reclutamiento
  BPR[i]=sum(N*wedad*Om*exp(-Z*dts)) #Biomasa desovante
  
  
}

B0=BPR[1]
alfa=4*h*R0/(5*h-1);
beta=(1-h)/(5*h-1)*B0;


BPR_eq=alfa*BPR-beta
Reclu=alfa*BPR_eq/(beta+BPR_eq)
YPR_eq=YPR*Reclu





par(mfrow = c(2, 2))

id=which(YPR_eq==max(YPR_eq))

plot(Fref,YPR,type="l",col="green",lwd=2)
lines(Fref,YPR_eq,col="red",lwd=2)
legend(0.6,250,c("h=1",paste("h=",h)),bty="n",col=c("green","red"),lty=1,lwd=2)
abline(v=Fref[id],lty=2)
abline(v=2*M,lty=2,col="red")

plot(Fref,BPR/BPR[1],type="l",col="green",lwd=2,ylim=c(0,1),ylab="B/B0")
lines(Fref,BPR_eq/BPR_eq[1],col="red",lwd=2)
legend(0.6,0.8,c("h=1",paste("h=",h)),bty="n",col=c("green","red"),lty=1,lwd=2)
abline(v=Fref[id],lty=2)
abline(v=2*M,lty=2,col="red")

plot(BPR/BPR[1],YPR,type="l",col="green",lwd=2,xlab="B/B0",ylab="Y")
lines(BPR_eq/BPR_eq[1],YPR_eq,col="red",lwd=2)
legend(0.6,250,c("h=1",paste("h=",h)),bty="n",col=c("green","red"),lty=1,lwd=2)
abline(v=BPR_eq[id]/BPR_eq[1],lty=2)


plot(edad,Sel,type="l",col="blue",ylab="Proporcion",main="Selectividad-Madurez",lwd=2)
lines(edad,Om,col="red",lwd=2)
legend(7,0.8,c("Selectividad","madurez"),bty="n",col=c("blue","red"),lty=1,lwd=2)


print(data.frame(F=Fref,BPR=BPR,YPR=YPR,BPR_eq=BPR_eq,YPR_eq=YPR_eq,pB0=BPR/BPR[1],pB0eq=BPR_eq/BPR_eq[1]))

#########################

plot(Fref,BPR/BPR[1],type="l",col="red",lwd=2,ylim=c(0,1),xlab='Mort por pesca',ylab="Rendimiento/Biomasa",main="Análisis por recluta")
lines(Fref,YPR/max(YPR),col="green",lwd=2)
legend("topright",c("B","Y"),col=c("red","green"),lty=1,bty="n",lwd=2)

plot(Fref,BPR_eq/BPR_eq[1],type="l",col="red",lwd=2,ylim=c(0,1),xlim=c(0,1),xlab='Mort por pesca',ylab="Rendimiento/Biomasa", main="Análisis de equilibrio por recluta")
lines(Fref,YPR_eq/max(YPR_eq),col="green",lwd=2)
abline(h=0.2,v=0.165,lty=2,col="red")
abline(h=0.4,v=0.295,lty=2,col="red")
abline(h=0.975,lty=2,col="green")
abline(h=0.931,lty=2,col="green")
arrows(x0 = 0.28,
       x1 = 0.18,
       y0 = 0.1,col="red",length = 0.1) #,y1 = 90) 
arrows(x0 = 0.5,
       y0 = 0.931,y1=0.975,col="green",length = 0.05) #,y1 = 90) 
legend("topright",c("B","Y"),col=c("red","green"),lty=1,bty="n",lwd=2)

plot(BPR,Reclu,type="l",col="black",lwd=2,ylim=c(0,1),ylab="%R0", main="Relación S/R")

plot(BPR_eq/BPR_eq[1],YPR_eq/max(YPR_eq),type="l",col="green",lwd=2,ylim=c(0,1),xlab='B/B0',ylab="Y/Ymax", main="Producción en equilibrio")
abline(v=0.334462701,lty=2,col="red")
text(0.3,0.05,paste("Brms/B0=0.33"),col="red",cex=1.1)
############################################
par(mfrow = c(1, 1))
plot(Fref,BPR_eq/BPR_eq[1],type="l",col="green",lwd=2,ylim=c(0,1),xlab='Mort por pesca',ylab='B/B0')

A50=4.22# edad completamente reclutada # 50%
A95=5
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50))) # selectividad
for (i in 1:length(Fref)){
  
  
  F=Fref[i]*Sel
  Z=F+M
  
  # dinámica para B0
  N[1]=R0
  
  for (a in 2:Tmax){
    N[a]=N[a-1]*exp(-M)
  }
  
  N[Tmax]=N[Tmax]/(1-exp(-M)) # número del grupo plus
  B0=sum(N*wedad*Om*exp(-M*dts))
  
  
  # dinámica en F
  
  N[1]=R0
  
  for (a in 2:Tmax){
    N[a]=N[a-1]*exp(-Z[a-1])
  }
  
  N[Tmax]=N[Tmax]/(1-exp(-Z[Tmax]))
  
  C=F/Z*N*(1-exp(-Z))
  
  YPR[i]=sum(C*wedad) #Reclutamiento
  BPR[i]=sum(N*wedad*Om*exp(-Z*dts)) #Biomasa desovante
  
  
}

B0=BPR[1]
alfa=4*h*R0/(5*h-1);
beta=(1-h)/(5*h-1)*B0;


BPR_eq=alfa*BPR-beta
Reclu=alfa*BPR_eq/(beta+BPR_eq)
YPR_eq=YPR*Reclu


lines(Fref,BPR_eq/BPR_eq[1],lwd="2",col="red")
abline(v=0.3,lty=2,col="black")
abline(h=0.4,lty=2,col="black")
arrows(x0 = 0.4,
       y0 = 0.15,y1 = 0.3,col="black",length = 0.1) #,y1 = 90) 
legend("topright",c("41","29.7"),col=c("red","green"),lty=1,bty="n",lwd=2,title='Talla de captura')