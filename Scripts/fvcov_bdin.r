fvcov_bdin=function(data,parini){


h=1e-3
n=length(Y)
B=rep(0,1,n)
CPUEpred=rep(0,1,n)
J=matrix(0,n,length(parini))
correl=matrix(0,length(parini),length(parini))



for (j in 1:length(parini))
{  

par0=parini
par0[j]=parini[j]+h
    
K=exp(par0[1])
r=exp(par0[2])
sigma=exp(par0[3])
rho=exp(par0[4])

#-----------------------------------

B[1]=K

for (t in 2:n)
{ 
  B[t]=B[t-1] + r/rho*B[t-1]*(1-(B[t-1]/K)^rho) - Y[t-1]
}

 id=which(CPUEdat>0)
 q=exp(mean(log(CPUEdat[id]/B[id])))
 CPUEpred=q*B

J[,j]=1/h*(-n*log(0.5/(sigma*sqrt(6.283185)))+(1/sigma*(log(CPUEpred[id])-log(CPUEdat[id])))^2)


}

vcov=solve(t(J)%*%J)

cv_par=sqrt(diag(vcov))
sd_par=cv_par*exp(parini)

for (i in 1:length(parini)){
  for (j in 1:length(parini)){
    
    correl[i,j]=vcov[i,j]/(cv_par[i]*cv_par[j])
    
  }}


out=list(vcov=vcov,correl=correl)


}

