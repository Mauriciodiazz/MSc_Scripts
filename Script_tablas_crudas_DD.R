#### datos obtenidos desde la pagina, entonces debo filtrarlos de manera que se cree un DF por especie y dejar las columnas que me interesan

gbif_crudo<-read.csv2("E:/Maestria_DD/spp. records_DD/specialists_DD/0098854-230224095556074.csv", sep="\t", stringsAsFactors = F, quote = "")
length(gbif_crudo$gbifID)
which(is.na(gbif_crudo$decimalLatitude))
which(gbif_crudo$gbifID==259530407)

gbif_crudo$decimalLatitude<-as.numeric(gbif_crudo$decimalLatitude)
gbif_crudo$decimalLongitude<-as.numeric(gbif_crudo$decimalLongitude)

str(gbif_crudo)

gbif_crudo[15569,]
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
str(gbif_cr.colfil)
which(is.na(gbif_cr.colfil$decimalLatitude))

#2. Partir la base de datos tablas por especie
spp_vec<-levels(factor(gbif_cr.colfil$species))
length(spp_vec)
lista<-list()

library('stringr')
for (x in 1:length(spp_vec)) {
  spp<-subset(gbif_cr.colfil, species==spp_vec[x])
  lista[[x]]<-spp
  names(lista)[x]<-spp_vec[x]
  lista[[x]]$locality<- str_replace_all(lista[[x]]$locality,'#','num ')
  write.table(lista[[x]], paste("E:/Maestria_DD/spp. records_DD/specialists_DD/spec. tablas crudas_DD/", spp_vec[x],".txt", sep = ""), sep = "\t", row.names = F, dec=".", quote = F)
}

a2<- read.table(paste("E:/Maestria_DD/spp. records_DD/specialists_DD/spec. tablas crudas_DD/", 
                        spp.list.vec[21], ".txt", sep=""),
                  header = T, sep = "\t", dec = ".", quote = "", 
                  blank.lines.skip=F)




length(lista[[31]]$gbifID)
str(lista[[31]])
length(a2$gbifID)

head(a2)
which(is.na(a2$decimalLatitude))

b<-a2[5566,]
b[2,]<-lista[[31]][5566,]
b[3,]<-a2[12809,]
b[4,]<-lista[[31]][12809,]
b[5,]<-a2[366140,]
b[6,]<-lista[[31]][366140,]
b[7,]<-a2[560681,]
b[8,]<-lista[[31]][560681,]


#row.names para que no me establezca la primera columna con los nombres de las filas
#quote porque GBIF usa caracteres vacíos dentro de columnas, entonces no quiero que las guarde como caracteres porque se modifican los nombres de las columnas 


a<-subset(gbif_cr.colfil, species== "Catharus dryas")
write.table(a, "E:/Maestria_DD/spp. records_DD/specialists_DD/cata.txt", sep = "\t", row.names = F, quote = F)

a.read<-read.table("E:/Maestria_DD/spp. records_DD/specialists_DD/spec. tablas crudas_DD/Aphelocoma woodhouseii.txt", header = T, sep = "\t", dec=".", blank.lines.skip = F, quote = "")

length(a$gbifID) == length(a.read$gbifID)


which(gbif_crudo$gbifID==793140858)
a<-gbif_crudo[which(gbif_crudo$gbifID==793140858),]
