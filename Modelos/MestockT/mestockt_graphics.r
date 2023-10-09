
rm(list=ls()) # erasure all objects
#system('./mestockt -ind datos10.dat -nox')  # for model running


source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')


data <-read.rep('for_R.rep')
attach(data)


target=0.4


barplot(Desembarques[1,]/1000~Yrs, cex.lab=1, xlab="Año",ylab="Toneladas (miles)", col = "gray",ylim = c(0,max(Desembarques[1,]/1000)*1.1), 
        main="Desembarques", cex.main = 1)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()

#Ajustes--------------------------------------------------------------------------

par(mfrow = c(2, 2))
sim=20
  
max=max(CPUE);
ubi=which(CPUE[1,]>0)
plot(Yrs[ubi],(CPUE[1,ubi]),ylab="Indice",xlab="Año",main="CPUE",ylim = c(min((CPUE[which(CPUE>0)])),max((CPUE[which(CPUE>0)]))),
     pch=sim,type="b",cex=1.5)
lines(Yrs[ubi],(CPUE[2,ubi]),col="red",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(Yrs,(Desembarques[1,]),ylab="Toneladas",xlab="Año",main="Desembarque",ylim = c(min((Desembarques)),max((Desembarques))),
     pch=sim,type="b",cex=1.5)
lines(Yrs,(Desembarques[2,]),col="red",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

max=max(Surveys[1,]);
ubi=which(Surveys[1,]>0)
plot(Yrs[ubi],(Surveys[1,ubi]),ylab="Indice",xlab="Año",main="Crucero",ylim = c(min((Surveys[which(Surveys>0)])),max((Surveys[which(Surveys>0)]))),
     pch=sim,type="b",cex=1.5)
lines(Yrs,(Surveys[2,]),col="red",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)





par(mfrow = c(2, 1))
sim=20
min=min(Lmed_flo[2:3,])
max=max(Lmed_flo[2:3,])
plot(Lmed_flo[1,],Lmed_flo[2,],ylab="Talla promedio",xlab="Año",main="Talla promedio flota",ylim = c(min,max),
     pch=sim,type="b",cex=1.5)
lines(Lmed_flo[1,],Lmed_flo[3,],lwd=2,col="red")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


min=min(Lmed_srv[2:3,])
max=max(Lmed_srv[2:3,])
plot(Lmed_srv[1,],Lmed_srv[2,],ylab="Talla promedio",xlab="Año",main="Talla promedio Crucero",ylim = c(min,max),
     pch=sim,type="b",cex=1.5)
lines(Lmed_srv[1,],Lmed_srv[3,],lwd=2,col="red")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

#Comps_tallas_marginal----------------------------------------------------------------------
#par(mfrow = c(1, 2))
dl=0.5*(Tallas[2]-Tallas[1])


plot(Tallas-dl, Frec_marg_flo[1,],col="gray",type="s",lwd=5,xlab="Talla",ylab="Proporcion",main="Flota")
lines(Tallas,Frec_marg_flo[2,],col="red",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(Tallas-dl, Frec_marg_srv[1,],col="gray",type="s",lwd=5,xlab="Talla",ylab="Proporcion",main="Cruceros")
lines(Tallas,Frec_marg_srv[2,],col="red",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


#Comps_tallas_f----------------------------------------------------------------------

n=length(Lmed_flo[1,])
filas=4 
cols=4
dl=0.5*(Tallas[2]-Tallas[1])

par(mfrow = c(filas, cols))
for (i in 1:n)
{
  plot(Tallas-dl,Frecs_capt_obs[i,],main=paste(Lmed_flo[1,i]),type="s",lwd=3, col="gray",cex.main = 1.2, ylab="", xlab="Talla",cex.lab=1.2,
               ylim=c(0,max(c(Frecs_capt_obs,Frecs_capt_pred))))
  lines(Tallas,Frecs_capt_pred[i,],col="red",lwd=2)
}


#Comps_tallas_s----------------------------------------------------------------------

n=length(Lmed_srv[1,])
filas=2
cols=round(n/filas)
if (cols==0)
{cols=1}

par(mfrow = c(filas, cols))

for (i in 1:n)
{

  plot(Tallas-dl,Frecs_srv_obs[i,],main=paste(c(Lmed_srv[1,i]),"c"),type="s",lwd=3, col="gray",cex.main = 1.2,cex.lab=1.2,
       ylab="", xlab="Talla",ylim=c(0,max(c(Frecs_srv_obs,Frecs_srv_pred))))
  lines(Tallas,Frecs_srv_pred[i,],col="red",lwd=2)

}


#Selectividad----------------------------------------------------------------------

par(mfrow = c(1, 2))

if(length(Sel_f)>length(L_edad))
{matplot(L_edad,t(Sel_f),type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",
         main="Selectividad Flota",col="gray",ylim=c(0,1))}

if (length(Sel_f)==length(L_edad))
{plot(L_edad,Sel_f,type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",
      main="Selectividad Flota",col="gray",ylim=c(0,1))}


lines(L_edad,Madurez_edad,col="red",lwd = 2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)



plot(L_edad,Sel_srv,lwd=2,type="l",ylim=c(0, 1),xlab="Talla",ylab="Proporcion",main="Selectividad Crucero",col="gray")
lines(L_edad,Madurez_edad,col="red",lwd = 2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


#Poblacionales----------------------------------------------------------------------

par(mfrow = c(2, 1))

sdR=Reclutamientos[3,]
li=Reclutamientos[2,]-1.96*sdR
ls=Reclutamientos[2,]+1.96*sdR
nyrs=length(Biomasa_desovante)
plot(Yrs,Reclutamientos[2,],ylim = c(0,max(Reclutamientos[2,]+1.96*sdR)*1.01),type="l",
     ylab="Reclutamientos",xlab="Año",pch = 16,cex=1,lwd=2,cex.main=2)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])
polygon(x,y,col="#DCDCDC",border="#DCDCDC")
lines(Yrs,Reclutamientos[2,],lwd="2")
lines(Yrs,Reclutamientos[1,],lwd="2",col="red")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(Yrs,Dev_log_R,type="l",xlab="Años",ylab="Desvio log_R",lwd=2,col="black",cex.main=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h=0,col="red",lwd = 2)

#Analisis por recluta------------

tmax=length(L_edad)


if(length(Sel_f)>length(L_edad))
{nsel=length(Sel_f[,1])
Selec=Sel_f[nsel,]}

if(length(Sel_f)==length(L_edad))
{Selec=Sel_f}

ypr_out<-por_recluta_r(tmax,Selec,Madurez_edad,Peso_edad,target,h,M,dts)
attach(ypr_out)
Ftar=Ftar[1]
BPRtar=BPRtar[1]
YPRtar=YPRtar[1]

par(mfrow = c(1, 1))
plot(Fcr,Y/max(Y),type="l", col="green", lwd=2, main="Analisis por recluta", xlab="Mortalidad por pesca", ylab="BPR, YPR relativos",
     cex.lab=1.1,cex.main=1.3,ylim = c(0,1))
lines(Fcr,B/max(B), col="red", lwd=2)
lines(Ftar, BPRtar/max(B),type="p", lwd=5)
text(Ftar, 0,Ftar,cex=1.1)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h = target, lty = 2,lwd=1)
abline(v = Ftar, lty = 2,lwd=1)


#Biomasa y F con IC-----------------------------------------------------------

par(mfrow = c(2, 1))

li=Biomasa_desovante[1,]-1.96*Biomasa_desovante[2,]
ls=Biomasa_desovante[1,]+1.96*Biomasa_desovante[2,]
nyrs=length(Biomasa_desovante)

plot(Yrs,Biomasa_desovante[1,],ylim = c(0,max(Biomasa_desovante[1,]+1.96*Biomasa_desovante[2,])*1.01),type="l",
     ylab="Biomasa",xlab="Año",pch = 16,cex=1,lwd=2,main="Biomasa", cex.lab=1.2,cex.main=1.2)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Biomasa_desovante[1,],lwd=2)
B0=mean(Biomasa_desovante[1,]/SPR[1,])
abline(h = target*B0, col = "red",lty = 2,lwd=2)



li=Mort_F[1,]-1.96*Mort_F[2,]
ls=Mort_F[1,]+1.96*Mort_F[2,]
nyrs=length(Mort_F)

plot(Yrs,Mort_F[1,],ylim = c(0,max(Mort_F[1,]+1.96*Mort_F[2,])*1.01),type="l",
     ylab="F",xlab="Año",pch = 16,cex=1,lwd=2,main="Mortalidad por pesca",cex.lab=1.2,cex.main=1.2)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Mort_F[1,],lwd=2)
abline(h = Ftar, col = "red",lty = 2,lwd=2)


#SPR------------------------------------------------------

par(mfrow = c(2, 1))
plot(Yrs,SPR[1,],type="l",xlab="Años",ylab="B/B0",ylim=c(0, max(SPR)),lwd=2,col="black",main="Potencial reproductivo (B/B0)", cex.lab=1.2,cex.main=1.2)
lines(Yrs,SPR[2,],col="green",lwd=2)
legend("topright",c("Largo plazo", "Dinñmico"),lty =c(1,1,1),col=c("black","green"),lwd=2,bty="n")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h = target, col = "red",lty = 2,lwd=2)

Bvirgen=Biomasa_desovante[1,]/SPR[2,]
plot(Yrs,Bvirgen,type="l",ylim=c(0, max(Bvirgen)),xlab="Años",ylab="Biomasa",lwd=2,col="green",main="Biomasa desovante", cex.lab=1.2,cex.main=1.2)
lines(Yrs,Biomasa_desovante[1,],col="black",lwd=2)
legend("topright",c("Virgen", "Explotada"),lty =c(1,1,1),col=c("green","black"),lwd=2,bty="n")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h = target*B0, col = "red",lty = 2,lwd=2)
abline(h = B0, col = "green",lty = 2,lwd=2)



#Kobe---------------------------------------------------------------------



par(mfrow = c(1, 1))
BRMS=B0*target
FRMS=Ftar
nyrs=length(Yrs)
  
plot(SPR[1,]/target,Mort_F[1,]/FRMS,pch = 16,ylab="F/Frms",xlab="B/Brms", xlim = c(0,max(SPR[1,]/target)), ylim = c(0,max(Mort_F[1,]/FRMS)*1.5), 
     type="o",col="black",lty="dashed",main=paste("B/Brms=",round(SPR[1,][length(Yrs)]/target,2),
                                                  " F/Frms=",round(Mort_F[1,][length(Yrs)]/FRMS,2)))

polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(0,0,1,1),col="green") #verde
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(1,1,1.5*max(Mort_F[1,]/FRMS),1.5*max(Mort_F[1,]/FRMS)),col="yellow1") #amarillo
polygon(c(0,1,1,0),c(1,1,1.5*max(Mort_F[1,]/FRMS),1.5*max(Mort_F[1,]/FRMS)),col="tomato1") #rojo

lines(SPR[1,]/target,Mort_F[1,]/FRMS,pch = 16, type="o",col="black",lty="dashed")
lines(SPR[1,nyrs]/target,Mort_F[1,][length(Yrs)]/FRMS,type="p",col="blue",pch = 16,cex=2)
text(SPR[1,]/target*.95,Mort_F[1,]/FRMS,paste(Yrs),cex=0.8)


X0=Biomasa_desovante[1,nyrs]/BRMS
Y0=Mort_F[1,nyrs]/FRMS
cvF=Mort_F[2,nyrs]/Mort_F[1,nyrs];
cvB=Biomasa_desovante[2,nyrs]/Biomasa_desovante[1,nyrs];

arrows(X0,Y0-1.96*cvF*Y0,X0,Y0+1.96*cvF*Y0,
       length = 0.15, code = 3, angle = 90, lwd=2, col="blue")

arrows(X0-1.96*cvB*X0,Y0,X0+1.96*cvB*X0,Y0,
       length = 0.15, code = 3, angle = 90,lwd=2, col="blue")

box()




#Genera excel--------------------------------------------------------------------------------------

R0=Reclutamientos[1,1]
Variables=data.frame(Año=Yrs,Biomasa=Biomasa_desovante[1,],R_R0=Reclutamientos[2,]/R0,
                     Fcr=Mort_F[1,], F_Fmrs=Mort_F[1,]/Ftar,B_Brms=Biomasa_desovante[1,]/BRMS,
                     B_B0=SPR[1,], SPRdin=SPR[2,])



write.csv(round(Variables,2), 'Var_Pobl.csv',  
          row.names = F,)

write.table(round(Variables,2), 'Var_Pobl.txt',  
          row.names = F,)



Variables2=data.frame(Año=Yrs,Desembarques=Desembarques[1,],CPUErel=CPUE[1,],
                     B.Acustica=Surveys[1,])


write.csv(Variables2, 'Data_usada.csv',  
          row.names = F,)




