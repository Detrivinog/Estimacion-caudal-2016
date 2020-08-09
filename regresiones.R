setwd("C:/Users/David Esteban/Universidad/Semestres/VII Semestre/Geomática Aplicada/Trabajo final")

datos<-read.csv2("datos_caudal.csv")
datos$E=datos$E/1000
datos$P=datos$P/1000

sub1<-subset(datos,NAME=="SC") # 1 La Chorrera 
sub2<-subset(datos,NAME=="SB") # 2 Santa Barbara
sub3<-subset(datos,NAME=="Chin") # 3 Chingaza


# Función de distribución de probabilidad Lognormal
fdp<-function(x,m,s){
  f=exp(-((log(x)-m)^2)/2*(s^2))/((x*s)*sqrt(2*pi))
}

summary(sub1)
min1=subset(sub1,Caud<= 1.453)
max1=subset(sub1,Caud>= 1.743)
sub1$min.u<-rep(mean(min1$Caud),length(sub1$ANO))
sub1$min.sd<-rep(sd(min1$Caud),length(sub1$ANO))
sub1$max.u<-rep(mean(max1$Caud),length(sub1$ANO))
sub1$max.sd<-rep(sd(max1$Caud),length(sub1$ANO))


summary(sub2)
min2=subset(sub2,Caud<=0.7157)
max2=subset(sub2,Caud>=0.8627)
sub2$min.u<-rep(mean(min2$Caud),length(sub2$ANO))
sub2$min.sd<-rep(sd(min2$Caud),length(sub2$ANO))
sub2$max.u<-rep(mean(max2$Caud),length(sub2$ANO))
sub2$max.sd<-rep(sd(max2$Caud),length(sub2$ANO))


summary(sub3)
min3=subset(sub3,Caud<= 6.651)
max3=subset(sub3,Caud>= 8.085)
sub3$min.u<-rep(mean(min3$Caud),length(sub3$ANO))
sub3$min.sd<-rep(sd(min3$Caud),length(sub3$ANO))
sub3$max.u<-rep(mean(max3$Caud),length(sub3$ANO))
sub3$max.sd<-rep(sd(max3$Caud),length(sub3$ANO))


mod1=lm(log(min.u)~log(A*(P-E)),sub1)
summary(mod1)
mod1s=lm(log(min.sd)~log(A*(P-E)),sub1)
summary(mod1s)


Tr=function(tempo){
  tr=vector()
  a=length(tempo)
  for(i in 1:(a-1)){
    tr[i]=tempo[i+1]-tempo[i]
  }
  c=mean(tr)
  tr[a]=c
  return(tr)
}

min1$t1<-Tr(min1$ANO)
min2$t2<-Tr(min2$ANO)
min3$t3<-Tr(min3$ANO)

max1$T1<-Tr(max1$ANO)
max2$T2<-Tr(max2$ANO)
max3$T3<-Tr(max3$ANO)
4.240799
min1$mean=rep(mean(min1$Caud),length(min1$ANO))
mod1=lm(mean~log(A*(P-E)), min1)
summary(mod1)

log(min1$mean)

sapply(min1, sd)
Qp=1.311944+fdp(2,4.375,4.240799)*0.107087
