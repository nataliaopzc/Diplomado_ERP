rm(list=ls())
inputs=read.csv("lbpa_dat_Almeja.csv",sep=";",encoding = "latin1")
par(mfrow = c(3, 2))
barplot(inputs$LF1~inputs$Length,ylim=c(0,max(inputs$LF1)), ylab = "Frecuencia", xlab = "",cex.lab=1.5)
barplot(inputs$LF2~inputs$Length,ylim=c(0,max(inputs$LF1)), ylab = "Frecuencia", xlab = "",cex.lab=1.5)
barplot(inputs$LF3~inputs$Length,ylim=c(0,max(inputs$LF1)), ylab = "Frecuencia", xlab = "",,cex.lab=1.5)
barplot(inputs$LF4~inputs$Length,ylim=c(0,max(inputs$LF1)), ylab = "Frecuencia", xlab = "Talla (mm)",cex.lab=1.5)
barplot(inputs$LF5~inputs$Length,ylim=c(0,max(inputs$LF1)), ylab = "Frecuencia", xlab = "Talla (mm)",cex.lab=1.5)

inputs2=read.csv("libro1.csv",sep=";",encoding = "latin1")
lmed=mean(inputs2$Frec)

