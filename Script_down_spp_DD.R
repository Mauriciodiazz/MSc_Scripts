###Script para limpiar datos MSc
### Mauricio Díaz - INECOL 2022

###Necesito comparar la lista de especies biologicas actualizada con la de la AOS para ver si hay alguna que no corresponda

AOS_63<- read.csv("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/NACC_list_species.csv")
head(AOS_63)

BSA<- read.table("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Biol_endemic spp_act.txt", sep='\t', dec=".", header=T)
head(BSA)

tabla<- data.frame(matrix(ncol=3, nrow = length(BSA$AOS.7th.63.suppl)))
names(tabla)<- c(names(BSA)[1],names(BSA)[2], "result")
tabla[,1]<- BSA[,1]
tabla[,2]<- BSA[,2]

for (x in 1:length(BSA$Biological.species)) {
  for (i in 1:length(AOS_63$species)) { #orden
  if(tabla[x,2] == AOS_63[i,9]) {tabla[x,3]<-"1"}
  }}

which(is.na(tabla$result))
tabla[which(is.na(tabla$result)),]

#esta parte es para especies endemicas
####Descargar registros GBIF####
##Opcion 1 - A través del paquete rgbif
spp<-read.table("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Biol_endemic spp_vector.txt", sep= "\t")
str(spp)

install.packages("rgbif")
library(rgbif)
?occ_search

spp[94,1]<-"Oporornis tolmiei" #esta mal escrito en el articulo de navarro

#con esta linea descargo los registros, -148 porque hay un ??? correspondiente a M. momotus 
#PROBLEMA: solo deja descargar máximo un millón de registros
#Solucion: restringir los registros a méxico
spp_list<- occ_search(scientificName = spp[-148,1], hasCoordinate= TRUE, limit= 5, geometry = "POLYGON((-92.62372 13.04972,-85.5032 19.55639,-87.0378 23.11665,-97.53442 27.59766,-101.34021 31.64899,-105.57569 33.18358,-109.3201 34.59541,-114.53771 34.34988,-118.89596 33.92019,-119.08011 30.23716,-110.11809 20.6613,-99.43731 14.40016,-92.62372 13.04972))",
                      fields = c("institutionCode", "catalogNumber", "collectionCode", "country", "stateProvince", "county", "locality", "year", "month", "decimalLatitude", "decimalLongitude", "elevation","basisOfRecord", "species", "family", "genus", "coordinateUncertaintyInMeters")) 


# get the DOIs for citing these data properly:
gbif_citation(gbif_data)
# note: if you need or prefer only one DOI for the entire dataset, download the dataset directly from www.gbif.org and then import the .csv to R. It is very important to properly cite the data sources! GBIF is not a source, just a repository for many people who put in very hard work to collect these data and make them available

#Cantidad total de registros en gbif por especie
recs<- data.frame(matrix(ncol=3, nrow = length(spp_list)))
names(recs)<- c("especies", "regs. GBIF")
str(recs)

for (x in 1:length(spp_list)) {
  recs[x,1]<-names(spp_list)[x]
  recs[x,2]<-spp_list[[x]][[1]][4]
}
recs

write.table(recs2, "C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/recs.txt", sep="\t")

#ordenar tabla y hacer un conteo de registros
rec_sort<-recs[order(recs$`regs. GBIF`, decreasing = FALSE), ]
rec_sort[1,3]<-rec_sort[1,2]

head(rec_sort)

#hacer una columna acumulada de las especies
for (i in 2:length(recs2$especies)) { #esta parte inicia el bucle y le dice a R que vamos a crear una variable llamada i que va a tomar valores desde i=2 hasta i=50.
  rec_sort[i,3] <- rec_sort[i,2]+rec_sort[i-1,3] #aquí, para cada valor de i, R aplica esa operación
} #esto termina el bucle.
head(rec_sort)


#darle un dataframe a cada especie y guardarlo en un txt 
for (x in 1:length(spp_list)) {
  a<- spp_list[[x]][[3]]
  write.table(a, paste("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/spp. tablas/",names(spp_list)[x],".txt", sep = ""), sep = "\t")
}
#### datos obtenidos desde la pagina, entonces debo filtrarlos de manera que se cree un DF por especie y dejar las columnas que me interesan

gbif_crudo<-read.csv2("E:/Maestria_DD/spp. records_DD/specialists/0087776-230224095556074/0087776-230224095556074.csv", sep="\t", stringsAsFactors = F, quote = "")
length(gbif_crudo$gbifID)
names(gbif_crudo)

#1. Filtrar los nombres de las columnas que me interesan
cols<-c("gbifID", "institutionCode", "catalogNumber", "collectionCode", "countryCode", "stateProvince", "locality", "year", "month", "decimalLatitude", "decimalLongitude", "elevation","basisOfRecord", "family", "genus", "species", "infraspecificEpithet", "scientificName", "coordinateUncertaintyInMeters", "coordinatePrecision")

col.ubi<- rep(NA, length(cols)) #este bucle establece la ubicación de esas columnas en la base de datos cruda
for (x in 1:length(cols)) {
  col.ubi[x]<-which(names(gbif_crudo)==cols[x])
}

gbif_cr.colfil<-gbif_crudo[,col.ubi] #selecciona esas columnas en un objeto nuevo
length(gbif_cr.colfil$gbifID)
names(gbif_cr.colfil)


#2. Partir la base de datos tablas por especie
spp_vec<-levels(factor(gbif_cr.colfil$species))
length(spp_vec)
lista<-list()

for (x in 1:length(spp_vec)) {
  spp<-subset(gbif_cr.colfil, species==spp_vec[x])
  lista[[x]]<-spp
  names(lista)[x]<-spp_vec[x]
  write.table(lista[[x]], paste("E:/Maestria_DD/spp. records_DD/specialists/spec. tablas crudas_DD/", spp_vec[x],".txt", sep = ""), sep = "\t", row.names = F, quote = F)
} 
#row.names para que no me establezca la primera columna con los nombres de las filas
#quote porque GBIF usa caracteres vacíos dentro de columnas, entonces no quiero que las guarde como caracteres porque se modifican los nombres de las columnas 

####Revisar datos cargados####
library(dplyr)
#Rupornis magnirostris
B.mag<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Tablas crudas_PC/R. magnirostris_0187395-220831081235567.csv", sep="\t")

R.mag<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Rupornis magnirostris.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)

length(R.mag$gbifID)
length(B.mag$gbifID)

levels(factor((R.mag$species)))
levels(factor((B.mag$species)))

levels(factor((R.mag$scientificName)))
levels(factor((B.mag$scientificName)))

#filtrar columnas de B. mag, es la descargada sola de GBIF
B.mag.colfi<-B.mag[,cols]
names(B.mag.colfi)
#write.table(B.mag.colfi, "D:/Maestria/Endemic spp rec/Tablas crudas/Buteo magnirostris.txt", sep = "\t", row.names=F)
RB.mag<-rbind(R.mag, B.mag.colfi)
levels(factor((RB.mag$species)))
levels(factor((RB.mag$scientificName)))

write.table(RB.mag, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Rupornis magnirostris_merge.txt", sep = "\t", row.names=F, quote=F)

#subespecies?
levels(factor(data$infraspecificEpithet))


#Aratinga holochlorus
A.hol<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Aratinga holochlora.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
P.hol<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Tablas crudas_PC/P. holochlorus_0187589-220831081235567.csv", sep="\t")
names(A.hol)

which(A.hol$catalogNumber=="OBS1190574316")

A.hol[93913,]
#son iguales?
length(A.hol$gbifID)
length(P.hol$gbifID)

levels(factor(A.hol$species))
levels(factor(P.hol$species))

levels(factor(A.hol$scientificName))
levels(factor(P.hol$scientificName))

#filtrar columnas de P. hol, es la descargada sola de GBIF
P.hol.colfi<-P.hol[,cols]
names(P.hol.colfi)
write.table(P.hol.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Psittacara holochlorus.txt", sep = "\t", row.names=F, quote = F)
AP.hol<-rbind(A.hol, P.hol.colfi)
levels(factor((AP.hol$species)))
levels(factor((AP.hol$scientificName)))

write.table(AP.hol, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Psittacara holochlorus_merge.txt", sep = "\t", row.names=F, quote = F)

#Momotus complex
M.coe<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Momotus coeruliceps.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
M.les<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Momotus lessonii.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
M.mom<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Momotus momota.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)

levels(factor(M.coe$species))
levels(factor(M.coe$scientificName))

levels(factor(M.les$species))
levels(factor(M.les$scientificName))

levels(factor(M.mom$species))
levels(factor(M.mom$scientificName))

M.sp<-rbind(M.coe, M.les, M.mom)
levels(factor(M.sp$species))

length(M.mom$gbifID)+length(M.coe$gbifID)+length(M.les$gbifID)
length(M.sp$gbifID)

write.table(M.sp, "D:/Maestria_DD//Endemic spp rec_DD/Tablas crudas_DD/Momotus sp_merge.txt", sep = "\t", row.names=F, quote=F)

#Colaptes rubiginosus
C.rub<-read.table("D:/Maestria_DD//Endemic spp rec_DD//Tablas crudas_DD/Colaptes rubiginosus.txt",  header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
P.rub<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Tablas crudas_PC/P. rubiginosus_0194196-220831081235567.csv", sep="\t")

levels(factor(C.rub$species))
levels(factor(P.rub$species))

levels(factor(C.rub$scientificName))
levels(factor(P.rub$scientificName))

#filtrar columnas de P.rub, es la descargada sola de GBIF
P.rub.colfi<-P.rub[,cols]
names(P.hol.colfi)
write.table(P.hol.colfi, "D:/Maestria_DD//Endemic spp rec_DD//Tablas crudas_DD/Piculus rubiginosus.txt", sep = "\t", row.names=F, quote = F)
PC.rub<-rbind(C.rub, P.rub.colfi)
levels(factor((PC.rub$species)))
levels(factor((PC.rub$scientificName)))

write.table(PC.rub, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Colaptes rubiginosus_merge.txt", sep = "\t", row.names=F, quote = F)

#Chlorospingus flavopectus
C.fla<-read.table("D:/Maestria_DD//Endemic spp rec_DD//Tablas crudas_DD/Chlorospingus flavopectus.txt",  header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
C.opt<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Tablas crudas_PC/C. ophthalmicus_0195098-220831081235567.csv", sep="\t")

levels(factor(C.fla$species))
levels(factor(C.opt$species))

levels(factor(C.fla$scientificName))
levels(factor(C.opt$scientificName))

levels(factor(C.fla$infraspecificEpithet))
levels(factor(C.opt$infraspecificEpithet))

#filtrar columnas de P.rub, es la descargada sola de GBIF
C.opt.colfi<-C.opt[,cols]
names(C.opt.colfi)
write.table(C.opt.colfi, "D:/Maestria_DD//Endemic spp rec_DD//Tablas crudas_DD/Chlorospingus ophthalmicus.txt", sep = "\t", row.names=F, quote = F)
C.fla.opt<-rbind(C.fla, C.opt.colfi)
levels(factor((C.fla.opt$species)))
levels(factor((C.fla.opt$scientificName)))

write.table(C.fla.opt, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Chlorospingus flavopectus_merge.txt", sep = "\t", row.names=F, quote = F)


#### Elimino registros duplicados####

#debo crear un bucle, que abra cada uno de los archivos dentro de la carpeta Tablas crudas_DD,luego que en cada uno busque duplicados y que me guarde la tabla sin duplicados en la carpeta spp. sin dup_DD

#1. vector de especies endemicas
spp.endemic<- read.table("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Endemic spp. biolcon.txt", sep = "\t")
spp.endemic.vec<-spp.endemic[,1]

#voy a crear una tabla que me muestre cuántos registros se borraron por especie
spp.filtrado<-data.frame(matrix(ncol=6, nrow=length(spp.endemic.vec)))
names(spp.filtrado)<- c("biosp_end", "raw", "sindup", "sinbuffer", "decimal", "spp")

for (x in 1:length(spp.endemic.vec)) {
  data<- read.table(paste("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/", spp.endemic.vec[x], ".txt", sep=""), 
                    header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
  spp.filtrado[x,1]<-spp.endemic.vec[x]
  spp.filtrado[x,2]<-length(data[,1])
  dups<- duplicated(data[, c("decimalLongitude", "decimalLatitude")])
  sindups <- data[!dups, ]
  spp.filtrado[x,3]<- length(sindups[,1])
 #write.table(sindups, paste("D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/", spp.endemic.vec[x],"_sd.txt", sep = ""), sep = "\t", row.names = F, quote = F)
 print(spp.filtrado[x,])
}

## Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 5 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
a<-read.table(paste("D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/", 
                    spp.endemic.vec[3], "_sd.txt", sep=""),
              header = T, sep = "\t", dec = ".", quote = "", 
              blank.lines.skip=F)
ata.clean[sapply(a$decimalLatitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
            sapply(data.clean$Latitude, decimalplaces) > 2, ]

for (x in 1:length(spp.endemic.vec)) {
  data.sd<- read.table(paste("D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/", 
                             spp.endemic.vec[x], "_sd.txt", sep=""),
                       header = T, sep = "\t", dec = ".", quote = "", 
                       blank.lines.skip=F)
  spp.filtrado[x,6]<-paste(spp.endemic.vec[x], "_sd_th.txt", sep="")
  data.clean<-clean_dup(data.sd, 
                        longitude= "decimalLongitude" , 
                        latitude = "decimalLatitude", 
                        threshold = 0.008333) #0.008333=1km
  spp.filtrado[x,4]<-length(data.clean[,1])
data.clean.dec<- data.clean[sapply(data.clean$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                              sapply(data.clean$decimalLatitude, decimalplaces) > 2, ]
  spp.filtrado[x,5]<-length(data.clean.dec[,1])
write.table(data.clean.dec, paste("D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/", spp.endemic.vec[x], "_sd_th.txt", sep = ""), sep = "\t", row.names = F, quote = F)
  print(spp.filtrado[x,])
}

####esta parte es para filtrar registros por cambios taxonómicos antes no detectados. Es decir, esta parte hace lo mismo de arriba pero por separado

####Aphelocoma coerulescens####
#esta especie fue descrita en 1944, pero luego se reconocieron algunos morfos como especies independientes, entre esos A. coerulescens, A. californica y A. woodhouseii; las dos primeras son especies endémicas de florida y california y la tercera es la especie que contiene los registros mexicanos (con las subespecies que componen la especie filogenética de Navarro). Los registros de GBIF descargados de A. coerulescens contienen los puntos de méxico, entonces voy a unir las dos tablas (A.coe+A.woo)

A.coe<-read.table("D:/Maestria_DD//Endemic spp rec_DD//Tablas crudas_DD/Aphelocoma coerulescens.txt",  header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)
A.woo<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Complementos_tablas_crudas_PC/A. woodhouseii_0216721-220831081235567.csv", sep="\t")

levels(factor(A.coe$species))
levels(factor(A.woo$species))

levels(factor(A.coe$scientificName))
levels(factor(A.woo$scientificName))

levels(factor(A.coe$infraspecificEpithet))
levels(factor(A.woo$infraspecificEpithet))

#filtrar columnas de A.woo, es la descargada sola de GBIF
A.woo.colfi<-A.woo[,cols]
names(A.woo.colfi)
write.table(A.woo.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Aphelocoma woodhouseii.txt", sep = "\t", row.names=F, quote = F)
A.coe.opt<-rbind(A.coe, A.woo.colfi)
levels(factor((A.coe.opt$species)))
levels(factor((A.coe.opt$scientificName)))

write.table(A.coe.opt, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Aphelocoma woodhouseii_merge.txt", sep = "\t", row.names=F, quote = F)

#Limpieza geografica
#1. duplicados
#data=A.coe.opt (el merge)

A.coe.dups<- duplicated(A.coe.opt[, c("decimalLongitude", "decimalLatitude")])
A.coe.sindups <- A.coe.opt[!A.coe.dups, ]

#write.table(A.coe.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Aphelocoma woodhouseii_merge_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=A.coe.sindups
A.coe.sindups$decimalLatitude<- as.numeric(A.coe.sindups$decimalLatitude)
A.coe.sindups$decimalLongitude<- as.numeric(A.coe.sindups$decimalLongitude)

data.clean.A.coe<-clean_dup(A.coe.sindups, 
                        longitude= "decimalLongitude", 
                        latitude = "decimalLatitude", 
                        threshold = 0.008333) #0.008333=1km

data.clean.dec<- data.clean.A.coe[sapply(data.clean.A.coe$Longitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                sapply(data.clean.A.coe$Latitude, decimalplaces) > 2, ]
  
#write.table(data.clean.dec, paste("D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/", spp.endemic.vec[x], "_sd_th.txt", sep = ""), sep = "\t", row.names = F, quote = F)

####Cynanthus doubledayi####
#esta especie ya es aceptada, por lo tanto está separada como un linaje biológico. lo que haré será filtrar solo los datos descargado de GBIF. Antes C. latirostris doubledayi

C.dou<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Complementos_tablas_crudas_PC/C. doubledayi_0240672-220831081235567.csv", sep="\t")

levels(factor(C.dou$species))
levels(factor(C.dou$scientificName))
levels(factor(C.dou$infraspecificEpithet))

#filtrar columnas de C.dou, es la descargada sola de GBIF
C.dou.colfi<-C.dou[,cols]
names(C.dou.colfi)
write.table(C.dou.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Cynanthus doubledayi.txt", sep = "\t", row.names=F, quote = F)

#Limpieza geografica
#1. duplicados
#data=C.dou.opt (el merge)

C.dou.dups<- duplicated(C.dou.colfi[, c("decimalLongitude", "decimalLatitude")])
C.dou.sindups <- C.dou.colfi[!C.dou.dups, ]

write.table(C.dou.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Cynanthus doubledayi_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=C.dou.sindups
C.dou.sindups$decimalLatitude<- as.numeric(C.dou.sindups$decimalLatitude)
C.dou.sindups$decimalLongitude<- as.numeric(C.dou.sindups$decimalLongitude)

data.clean.C.dou<-clean_dup(C.dou.sindups, 
                            longitude= "decimalLongitude", 
                            latitude = "decimalLatitude", 
                            threshold = 0.008333) #0.008333=1km

data.clean.dec<- data.clean.C.dou[sapply(data.clean.C.dou$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                    sapply(data.clean.C.dou$decimalLatitude, decimalplaces) > 2, ]

write.table(data.clean.dec, "D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/Cynanthus doubledayi_sd_th.txt", sep = "\t", row.names = F, quote = F)


####Euphonia godmani####
#esta especie ya es aceptada, por lo tanto está separada como un linaje biológico. lo que haré será filtrar solo los datos descargado de GBIF. Antes Euphonia affinis godmani

E.god<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Complementos_tablas_crudas_PC/E. godmani_0240787-220831081235567.csv", sep="\t")

levels(factor(E.god$species))
levels(factor(E.god$scientificName))
levels(factor(E.god$infraspecificEpithet))

#filtrar columnas de A.woo, es la descargada sola de GBIF
E.god.colfi<-E.god[,cols]
names(E.god.colfi)
write.table(E.god.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Euphonia godmani.txt", sep = "\t", row.names=F, quote = F)

#Limpieza geografica
#1. duplicados
#data=E.god.opt (el merge)

E.god.dups<- duplicated(E.god.colfi[, c("decimalLongitude", "decimalLatitude")])
E.god.sindups <- E.god.colfi[!E.god.dups, ]

write.table(E.god.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Euphonia godmani_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=E.god.sindups
E.god.sindups$decimalLatitude<- as.numeric(E.god.sindups$decimalLatitude)
E.god.sindups$decimalLongitude<- as.numeric(E.god.sindups$decimalLongitude)

data.clean.E.god<-clean_dup(E.god.sindups, 
                            longitude= "decimalLongitude", 
                            latitude = "decimalLatitude", 
                            threshold = 0.008333) #0.008333=1km

data.clean.dec<- data.clean.E.god[sapply(data.clean.E.god$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                    sapply(data.clean.E.god$decimalLatitude, decimalplaces) > 2, ]

write.table(data.clean.dec, "D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/Euphonia godmani_sd_th.txt", sep = "\t", row.names = F, quote = F)

####Junco bairdi####
#esta especie ya es aceptada, por lo tanto está separada como un linaje biológico. lo que haré será filtrar solo los datos descargado de GBIF. Antes Euphonia affinis godmani

J.bai<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Complementos_tablas_crudas_PC/J. bairdi_0250689-220831081235567.csv", sep="\t")

levels(factor(J.bai$species))
levels(factor(J.bai$scientificName))
levels(factor(J.bai$infraspecificEpithet))

#filtrar columnas de J.bai, es la descargada sola de GBIF
J.bai.colfi<-J.bai[,cols]
names(J.bai.colfi)
write.table(J.bai.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Junco bairdi.txt", sep = "\t", row.names=F, quote = F) #RECORDAR CAMBIAR LOS NOMBRES!!!!

#Limpieza geografica
#1. duplicados
#data=J.bai.opt (el merge)

J.bai.dups<- duplicated(J.bai.colfi[, c("decimalLongitude", "decimalLatitude")])
J.bai.sindups <- J.bai.colfi[!J.bai.dups, ]

write.table(J.bai.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Junco bairdi_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=J.bai.sindups
J.bai.sindups$decimalLatitude<- as.numeric(J.bai.sindups$decimalLatitude)
J.bai.sindups$decimalLongitude<- as.numeric(J.bai.sindups$decimalLongitude)

data.clean.J.bai<-clean_dup(J.bai.sindups, 
                            longitude= "decimalLongitude", 
                            latitude = "decimalLatitude", 
                            threshold = 0.008333) #0.008333=1km

data.clean.dec<- data.clean.J.bai[sapply(data.clean.J.bai$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                    sapply(data.clean.J.bai$decimalLatitude, decimalplaces) > 2, ]

write.table(data.clean.dec, "D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/Junco bairdi_sd_th.txt", sep = "\t", row.names = F, quote = F)

####Momotus coeruliceps####
#esta especie ya es aceptada, por lo tanto está separada como un linaje biológico. lo que haré será filtrar solo los datos descargado de GBIF.

M.coe<-read.table("D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Momotus coeruliceps.txt", header = T, sep = "\t", dec = ".", quote = "", blank.lines.skip=F)

levels(factor(M.coe$species))
levels(factor(M.coe$scientificName))
levels(factor(M.coe$infraspecificEpithet))

#filtrar columnas de M.coe, es la descargada sola de GBIF
M.coe.colfi<-M.coe[,cols]
names(M.coe.colfi)
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(M.coe.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Momotus coeruliceps.txt", sep = "\t", row.names=F, quote = F) 

#Limpieza geografica
#1. duplicados
#data=M.coe.opt (el merge)

M.coe.dups<- duplicated(M.coe.colfi[, c("decimalLongitude", "decimalLatitude")])
M.coe.sindups <- M.coe.colfi[!M.coe.dups, ]
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(M.coe.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Momotus coeruliceps_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=M.coe.sindups
M.coe.sindups$decimalLatitude<- as.numeric(M.coe.sindups$decimalLatitude)
M.coe.sindups$decimalLongitude<- as.numeric(M.coe.sindups$decimalLongitude)

data.clean.M.coe<-clean_dup(M.coe.sindups, 
                            longitude= "decimalLongitude", 
                            latitude = "decimalLatitude", 
                            threshold = 0.008333) #0.008333=1km

data.clean.M.coe.dec<- data.clean.M.coe[sapply(data.clean.M.coe$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                    sapply(data.clean.M.coe$decimalLatitude, decimalplaces) > 2, ]
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(data.clean.M.coe.dec, "D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/Momotus coeruliceps_sd_th.txt", sep = "\t", row.names = F, quote = F)

####Phaethornis mexicanus####
#esta especie ya es aceptada, por lo tanto está separada como un linaje biológico. lo que haré será filtrar solo los datos descargado de GBIF. Antes Phaethornis longirostris mexicanus y P. l. griseoventer

P.mex<-read.csv2("C:/Users/Mauricio Diaz/Documents/Maestria/Endemic spp rec/Complementos_tablas_crudas_PC/P. mexicanus_0251976-220831081235567.csv", sep="\t")

levels(factor(P.mex$species))
levels(factor(P.mex$scientificName))
levels(factor(P.mex$infraspecificEpithet))

#filtrar columnas de P.mex, es la descargada sola de GBIF
P.mex.colfi<-P.mex[,cols]
names(P.mex.colfi)
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(P.mex.colfi, "D:/Maestria_DD/Endemic spp rec_DD/Tablas crudas_DD/Phaethornis mexicanus.txt", sep = "\t", row.names=F, quote = F) 

#Limpieza geografica
#1. duplicados
#data=P.mex.opt (el merge)

P.mex.dups<- duplicated(P.mex.colfi[, c("decimalLongitude", "decimalLatitude")])
P.mex.sindups <- P.mex.colfi[!P.mex.dups, ]
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(P.mex.sindups, "D:/Maestria_DD/Endemic spp rec_DD/spp. sin dup_DD/Phaethornis mexicanus_sd.txt", sep = "\t", row.names = F, quote = F)

## 2. Eliminar valores muy cercanos y con menos de 2 decimales en lon y lat
## Quitar los registros bajo un umbral o distancia espacial. Por ejemplo 1 km

library(hsi)
resol <- 0.008333 # 1 km

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
# data.sd=P.mex.sindups
P.mex.sindups$decimalLatitude<- as.numeric(P.mex.sindups$decimalLatitude)
P.mex.sindups$decimalLongitude<- as.numeric(P.mex.sindups$decimalLongitude)

data.clean.P.mex<-clean_dup(P.mex.sindups, 
                            longitude= "decimalLongitude", 
                            latitude = "decimalLatitude", 
                            threshold = 0.008333) #0.008333=1km

data.clean.P.mex.dec<- data.clean.P.mex[sapply(data.clean.P.mex$decimalLongitude, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                                          sapply(data.clean.P.mex$decimalLatitude, decimalplaces) > 2, ]
#RECORDAR CAMBIAR LOS NOMBRES!!!!
write.table(data.clean.P.mex.dec, "D:/Maestria_DD/Endemic spp rec_DD/spp. clean_DD/Phaethornis mexicanus_sd_th.txt", sep = "\t", row.names = F, quote = F)

##Final del script!####

