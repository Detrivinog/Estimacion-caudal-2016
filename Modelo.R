setwd("C:/Users/David Esteban/Google Drive/Materias cursadas/Geomática Aplicada/Trabajo final")

#Paquetes utilizados 

library(sp)
library(gstat)
library(rgdal)
require(automap)
library(maptools)
library(rgeos)
library(plyr)

#Datos

series<-read.csv2("series.csv")
preci<-subset(series,DATA=="PRECIP")
temp=subset(series,DATA=="TEMP")

for(i in 1:dim.data.frame(temp)[1]){
  temp$I[i]=sum((temp[i,3:14]^1.514)/(5^1.514))
}
head(temp)
#Thornwite en mm/mes
EvP<-function(t,i){
  I=i
  a=0.49+1.79e-2*I-7.71e-5*(I^2)+6.75e-7*(I^3)
  ep=16*(10*t/I)^a
  return(ep)
}
dim.data.frame(temp)
EP=data.frame()

for(i in 1:dim.data.frame(temp)[1]){
  for(j in 1:dim.data.frame(temp)[2]){
    if(j==1){
      EP[i,j]=temp$ANO[i]
    }
    if(j==2){
      EP[i,j]=temp$DEP[i]
    }
   if(j>2 & j<=14){
     EP[i,j]=EvP(temp[i,j],temp$I[i])
   }
    if(j==15){
      EP[i,j]=sum(EP[i,3:14])
    }
  if(j>15){
    EP[i,j]=temp[i,j]
   }
  }
}
colnames(EP)=colnames(temp)
head(EP)   
#Evapotranspiración potencial 
write.csv2(EP,"EP1.csv")
#Shapes de las cuencas 

ching<-readOGR(".","CHING")
sc<-readOGR(".","SC")
sb<-readOGR(".","SB")

plot(ching, axes=T, col="gray74")
plot(sc)

gArea(ching)/10000
#Georreferenciación de las bases de datos
coordinates(EP)<-~X+Y
proj4string(EP)<-CRS("+init=epsg:4326")
EP<- spTransform(EP, crs.new)

coordinates(preci)<-~X+Y
proj4string(preci)<-CRS("+init=epsg:4326")
precip<- spTransform(preci, crs.new)

plot(sub1,axes=TRUE, main="Cuencas", cex.axis=.95)
plot(sc, col="blue", add=T)
plot(sb,col="green", add=T)
plot(ching,col="red",add=T)

sub1=subset(precip)
plot(sub1)
sub1
krig1= autoKrige(ANUAL~1, sub1,sc)
krig1$krige_output@data$var1.pred
p1<-plot(krig1)

sub1@data[1]
names(sub1)[1:14]
as.formula(paste("ENE","~","1"))

#Función de interpolación, dado los datos, el shape de salida, la formula con la que se maneja el kriging, 
# y variable del departamento. 
interpol<-function(data, shape,tempo, Departamento){
  sub1=subset(data, DEP==Departamento)
  t=length(tempo)
  interpol=data.frame()
  for(i in 1:t){
    sub=subset(sub1,ANO==c(0,tempo[i]))
    nombre=names(sub)
    for(j in 1:14){
      if(j<3){
        interpol[i,j]=sub@data[i,j]
      }else{
        form=as.formula(paste(nombre[j],"~","1"))
       krig= autoKrige(form, sub, shape)
       interpol[i,j]=krig$krige_output@data$var1.pred
      }
    }
  }
  colnames(interpol)=names(sub1)[1:14]
  return(interpol)
}

    ############################
    ## Ecuación Budyko (1974) ##
    ############################

Eva=function(ep,p){
  e=sqrt(((ep*p)*tanh(p/ep))*(1-cosh(ep/p)+sinh(ep/p)))
}


  #############################
  #####  San Carlos (1)  ######
  #############################


tempo=seq(1980,2015, 1)
scEP<-interpol(EP,sc,tempo,"ANT")
scP<-interpol(precip,sc,tempo,"ANT")

head(scEP)
scA<-gArea(sc)
scE<-Eva(scEP[3:14],scP[3:14])
head(scE)
for (i in 1:36) {
  scP$P[i]=sum(scP[i,3:14])
}
for (i in 1:36) {
  scE$E[i]=sum(scE[i,1:12])
}
#Caudal de la cuanca en m3*seg
scCaudal<-((scP[3:14]-scE[1:12])/1000)*scA/(30*24*60*60)
head(scCaudal)
scCaudal$name=rep("SC",length(tempo))
scCaudal$A<-rep(scA,length(tempo))
scCaudal$ANO=tempo
  ###############################
  #####  Santa Barbara (2)  #####
  ###############################

tempo=seq(1980,2015, 1)
sbEP<-interpol(EP,sb,tempo,"ANT")
sbP<-interpol(precip,sb,tempo,"ANT")
for (i in 1:36) {
  sbP$P[i]=sum(sbP[i,3:14])
}
for (i in 1:36) {
  sbE$E[i]=sum(sbE[i,1:12])
}

sbE<-Eva(sbEP[3:14],sbP[3:14])

sbA<-gArea(sb)
#Caudal de la cuanca en m3*seg
sbCaudal<-((sbP[3:14]-sbE[1:12])/1000)*sbA/(30*24*60*60)
head(sbCaudal)
sbCaudal$name2=rep("SB",length(tempo))
sbCaudal$A2<-rep(sbA,length(tempo))
sbCaudal$ANO=tempo
head(scCaudal)
  ############################
  #####  Chingaza  (3)  ######
  ############################

tempo=seq(1980,2015, 1)
chinEP<-interpol(EP,ching,tempo,"CUN")
chinP<-interpol(precip,ching,tempo,"CUN")
for (i in 1:36) {
  chinP$P[i]=sum(chinP[i,3:14])
}

chinE<-Eva(chinEP[3:14],chinP[3:14])
for (i in 1:36) {
  chinE$E[i]=sum(chinE[i,1:12])
}
chinA<-gArea(ching)
#Caudal de la cuanca en m3*seg
chinCaudal<-((chinP[3:14]-chinE[1:12])/1000)*chinA/(30*24*60*60)
chinCaudal$name3=rep("Chin",length(tempo))
chinCaudal$A3=rep(chinA,length(tempo))
chinCaudal$ANO=tempo
head(chinCaudal)


extremos=function(caudal){
  minimo=vector()
  sd.min=vector()
  maximo=vector()
  sd.max=vector()
  medio=vector()
  sd.medio=vector()
  for (i in 1:36) {
    x=as.vector(as.matrix(caudal[i,1:12]))
    q<-quantile(x)
    minimo[i]=mean(subset(x,x<q[[2]]&x>0))
    sd.min[i]=sd(subset(x,x<q[[2]]))
    maximo[i]=mean(subset(x,x>q[[4]]))
    sd.max[i]=sd(subset(x,x>q[[4]]))
    medio[i]=mean(x)
    sd.medio[i]=sd(x)
  }
  extrem=data.frame(medio,sd.medio,minimo,sd.min,maximo,sd.max)
  return(extrem)
}


c3=cbind(chinCaudal,extremos(chinCaudal))
c3$E=chinE$E
c3$P=chinP$P
c2=cbind(sbCaudal,extremos(sbCaudal))
c2$E=sbE$E
c2$P=sbP$P
c1=cbind(scCaudal,extremos(scCaudal))
c1$E=scE$E
c1$P=scP$P

mean(c1$medio)
mean(c1$maximo)
mean(c1$minimo)

mean(c2$medio)
mean(c2$maximo)
mean(c2$minimo)

mean(c3$medio)
mean(c3$maximo)
mean(c3$minimo)

mean(c1$E)
mean(c2$E)
mean(c3$E)

mean(c1$P)
mean(c2$P)
mean(c3$P)


mod1=lm(log(maximo)~log(medio),c1)
summary(mod1)
mod2=lm(log(sd.max)~log(medio),c1)
summary(mod2)
mod3=lm(log(minimo)~log(A*(P-E)),c1)
summary(mod3)
c1$sd.min
mod4=lm(log(sd.min)~log(A*(P-E)),c1)
summary(mod4)

mod1=lm(log(maximo)~log(medio),c2)
summary(mod1)
mod2=lm(log(sd.max)~log(medio),c2)
summary(mod2)
mod3=lm(log(minimo)~log(medio),c2)
summary(mod3)
c1$sd.min
mod4=lm(log(sd.min)~log(medio),c2)
summary(mod4)

mod1=lm(log(maximo)~log(medio),c3)
summary(mod1)
mod2=lm(log(sd.max)~log(medio),c3)
summary(mod2)
mod3=lm(log(minimo)~log(medio),c3)
summary(mod3)
c1$sd.min
mod4=lm(log(sd.min)~log(medio),c3)
summary(mod4)
c2$E==c3$E

