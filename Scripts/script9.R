rm(list=ls()) # comando para limpiar


# Parametros biologicos
M=0.2
Tmax=12
A50=5
A95=7
Fref=2.5*M
R0=1
dtf=0.5
dts=0.8
Loo=150
k=0.15
t0=-0.5
aw=0.001
bw=3.0
A50m=6
A95m=7


# vectores biologicos
edad=seq(1,Tmax)
Ledad=Loo*(1-exp(-k*(edad-t0)))
wedad=aw*Ledad^bw
N=rep(0,Tmax)
Om=1/(1+exp(-log(19)*(edad-A50m)/(A95m-A50m)))
Sel=1/(1+exp(-log(19)*(edad-A50)/(A95-A50)))
F=Fref*Sel
Z=F+M

# dinámica
N[1]=R0

for (a in 2:Tmax){
  N[a]=N[a-1]*exp(-Z[a-1])
}

N[Tmax]=N[Tmax]/(1-exp(-Z[Tmax]))

C=F/Z*N*(1-exp(-Z))

par(mfrow = c(1, 2))

plot(edad,Sel,type="l",lwd=2)
lines(edad,Om,lwd=2,col="green")
abline(h=0.5,lty=2,col="red")

Y=sum(C*wedad)
BD=sum(N*wedad*Om*exp(-Z*dts))

barplot(N~edad)
barplot(C~edad, add = T,col = "red")

print(paste("YPR=", round(Y,2), "BPR=", round(BD,2)))









