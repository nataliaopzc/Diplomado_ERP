

# Parametros biologicos
M=0.33
Tmax=12
A50=2.4
A95=3
Fref=0 #para calculat el ratio de agotamaiento (Biomasa actual /bomasa virginial (F=0)). Ahí te da el % de biomasa en el agua y puedes ver qué F te sirve

#Fref=2.5*M #para calculat el ratio de agotamaiento (Biomasa actual /bomasa virginial (F=0)). Ahí te da el % de biomasa en el agua y puedes ver qué F te sirve
R0=1
dtf=0.5
dts=0.5833
Loo=800.4
k=0.14
t0=-0.918
aw=0.001
bw=3.0
A50m=2.8
A95m=3.5


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

Y=sum(C*wedad) #Rendimiento
BD=sum(N*wedad*Om*exp(-Z*dts)) # Biomasa desovante

barplot(N~edad)
barplot(C~edad, add = T,col = "red")

print(paste("YPR=", round(Y,2), "BPR=", round(BD,2)))









