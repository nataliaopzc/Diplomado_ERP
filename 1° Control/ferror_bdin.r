ferror_bdin=function(data,parini){

K=exp(parini[1])
r=exp(parini[2])
sigma=exp(parini[3])

if(opt_p[1]==1){
  rho=exp(parini[4])}
  else{
    rho=p[1]}


n=length(Y)
B=rep(0,1,n)
CPUEpred=rep(0,1,n)

B[1]=K

for (t in 2:n)
{ 
  B[t]=B[t-1] + r/rho*B[t-1]*(1-(B[t-1]/K)^rho) - Y[t-1]
}

 id=which(CPUEdat>0)
 q=exp(mean(log(CPUEdat[id]/B[id])))
 CPUEpred=q*B

suma=sum((1/sigma*(log(CPUEpred[id])-log(CPUEdat[id])))^2)
fun=-n*log(0.5/(sigma*sqrt(6.283185)))+suma

out=suma
return(fun)

}

