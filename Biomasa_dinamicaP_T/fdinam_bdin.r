fdinam_bdin=function(data,parini){

K=exp(parini[1])
r=exp(parini[2])
sigma=exp(parini[3])

n=length(Y)
rho=p[1]
Biom=rep(0,1,n)
CPUEpred=rep(0,1,n)

Biom[1]=K

for (t in 2:n)
{ 
  Biom[t]=max(c(Biom[t-1] + r/rho*Biom[t-1]*(1-(Biom[t-1]/K)^rho) - Y[t-1],0.1))
}

 id=which(CPUEdat>0)
 q=exp(mean(log(CPUEdat[id]/Biom[id])))
 CPUEpred=q*Biom
 G=r/rho*Biom*(1-(Biom/K)^rho)
 Fmort=Y/Biom
 
salidas=data.frame(Y,CPUEdat,CPUEpred,Biom,G,Fmort) # salidas


return(salidas)

}

