rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=0.2
Tmax=12
A50=4 #5
A95=6 #7
R0=1
dts=0.8 #dt desove
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
Om=1/(1+exp(-log(19)*(edad-A50m)/(A95m-A50m)))
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50)))


Fref=seq(0,5*M,0.01) # genera un vector de valores de F
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

N[Tmax]=N[Tmax]/(1-exp(-M))
B0=sum(N*wedad*Om*exp(-M*dts))


# dinámica en F

N[1]=R0

for (a in 2:Tmax){
  N[a]=N[a-1]*exp(-Z[a-1])
}

N[Tmax]=N[Tmax]/(1-exp(-Z[Tmax]))

C=F/Z*N*(1-exp(-Z))

YPR[i]=sum(C*wedad)
BPR[i]=sum(N*wedad*Om*exp(-Z*dts))


}


B0=BPR[1]
alfa=4*h*R0/(5*h-1);
beta=(1-h)/(5*h-1)*B0;


BPR_eq=alfa*BPR-beta
Reclu=alfa*BPR_eq/(beta+BPR_eq)
YPR_eq=YPR*Reclu


par(mfrow = c(2, 2))


id=which(YPR_eq-max(YPR_eq)==0)

plot(Fref,YPR/max(YPR),type="l",col="green",lwd=2,xlab="Mortalidad por pesca",ylab="Rendimiento, Biomasa",
     main="Analisis por recluta (h=1)")
lines(Fref,BPR/BPR[1],type="l",col="red",lwd=2,ylim=c(0,1))
abline(v=Fref[id],lty=2)


plot(Fref,YPR_eq/max(YPR_eq),type="l",col="green",lwd=2,xlab="Mortalidad por pesca",ylab="Rendimiento, Biomasa",
     ylim=c(0,1),main=paste("Analisis de equilibrio (h=",h,")"))
lines(Fref,BPR_eq/BPR_eq[1],col="red",lwd=2)
abline(v=Fref[id],lty=2)
text(Fref[id]*1.1,0.1,round(Fref[id],2),col="red")
legend("topright",paste(c("Rendimiento","Biomasa")),col=c("green","red"),lty=1,bty="n")


plot(BPR_eq/BPR_eq[1],YPR_eq/max(YPR_eq),col="blue",lwd=2,type="l",xlab="Biomasa/B0",ylab="Rendimiento relativo",
     ylim=c(0,1),xlim=c(0,1),main=paste("Curva de producción Bmsy/B0=",round(BPR_eq[id]/BPR_eq[1],2)))
abline(v=BPR_eq[id]/BPR_eq[1],lty=2)


plot(edad,Sel,col="blue",lwd=2,type="l",xlab="Edad",main="Selectividad,madurez",ylab="Proporción")
lines(edad,Om,col="green",lwd=2)
legend("topleft",paste(c("Selectividad","Madurez")),col=c("blue","green"),lty=1,bty="n")


print(data.frame(F=Fref,BPR=BPR,YPR=YPR,BPR_eq=BPR_eq,YPR_eq=YPR_eq,pB0=BPR/BPR[1],pB0eq=BPR_eq/BPR_eq[1]))



