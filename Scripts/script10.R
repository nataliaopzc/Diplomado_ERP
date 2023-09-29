rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=0.2
Tmax=12
A50=5 
A95=6
R0=1
dts=0.8
Loo=150
k=0.15
t0=-0.5
aw=0.001
bw=3.0
A50m=6
A95m=7
h=0.75


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

N[Tmax]=N[Tmax]/(1-exp(-M)) # número del grupo pl
B0=sum(N*wedad*Om*exp(-M*dts))


# dinámica en F

N[1]=R0

for (a in 2:Tmax){
  N[a]=N[a-1]*exp(-Z[a-1])
}

N[Tmax]=N[Tmax]/(1-exp(-Z[Tmax]))

C=F/Z*N*(1-exp(-Z))

YPR[i]=sum(C*wedad) #Rendimiento
BPR[i]=sum(N*wedad*Om*exp(-Z*dts)) #Biomasa desovante


}

B0=BPR[1]
alfa=4*h*R0/(5*h-1);
beta=(1-h)/(5*h-1)*B0;


BPR_eq=alfa*BPR-beta #Biomasa en equilibrio a largo plazo
Reclu=alfa*BPR_eq/(beta+BPR_eq) #Reclutamiento
YPR_eq=YPR*Reclu #Rendiento en equilibrio a largo plazo





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



