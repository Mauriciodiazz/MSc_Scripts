###SIG - variables ambientales MSc
#Este script tiene los analisis de correlación
library(sf)
library(terra)
library(dplyr)

dem_ame<-rast("./WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_elev_America.tif")
ext(dem_ame) #este es el extent que contiene a méxico


# Preparacion de variables ------------------------------------------------


#Cargar capas mundiales y cortarlas al mismo extent


#####################################   slope ###
slope<-rast("./WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_slope.tif")
slope_ame<-crop(x=slope, y=dem_ame)
#writeRaster(slope_ame, "./WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_slope_Ame.tif", overwrite=TRUE)
#Slope no va a ir en los modelos!!!!


#####################################   Compound topographic index (CTI)##
#Esta variable esta hecha con la pendiente, hice una regresión entre ambas y me arrojó un R2 de 0.48... Entonces no la voy a usar
cti<-rast("./CTI/na_cti_9.tif")
cti
cti_1km<-resample(cti, dem_ame, method="bilinear") #esto cambia tambien el extent!!!!! Método bilinear porque son valores continuos

plot(cti_1km)
 ext(cti)
crs(cti_1km)==crs(dem_ame)
ext(cti_1km)==ext(dem_ame)
dim(cti_1km)==dim(dem_ame)

writeRaster(cti_1km, "./CTI/na_cti_9_Ame_1km.tif", overwrite=TRUE, names="cti_9_Ame_1km")

#cargar variables bioclimáticas
#con este bucle corto los raster
for (x in 1:19) {
  biox<-rast(paste("./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_",x,".tif", sep=""))
  biox_ame<-crop(x=biox, y=dem_ame)
  writeRaster(biox_ame, paste("./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/wc2.1_30s_bio_",x,"_Ame.tif", sep=""), overwrite=TRUE)
}


#####################################   Evapotranspiración###
#mes mas alto y mes mas bajo
#con este bucle corto los raster
for (x in 1:12) {
  evapx<-rast(paste("./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_",x,".tif", sep=""))
  evapx_ame<-crop(x=evapx, y=dem_ame)
  writeRaster(evapx_ame, paste("./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_Ame/wc2.1_30s_vapr_",x,"_Ame.tif", sep=""), overwrite=TRUE)
}

evap.stack<-rast(c(list.files(path = "./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_Ame/", pattern = "_Ame.tif$", full.names = TRUE)))
summary(evap.stack)

#con esta función obtengo estadísticos de los rasters
evap.stack %>%
  global(range, na.rm=TRUE)%>% #"min", "max", "sum", "prod", "range"
  arrange(desc(X2))%>%
  sum(X2,-X1)

arrange(a, .by_group = F)

rast("./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_1_MX.tif")
hist(as.data.frame(evap))

###como estas variables estan mes a mes, entonces debo unir todas las bbdd en un stack y de esta manera entonces calcular un nuevo raster con la media y otro con la mediana para luego escoger
evapx_ame_stack<- rast(list.files(path = "./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_Ame/", full.names = TRUE))

evapx_ame_mean<-app(evapx_ame_stack, mean)
plot(evapx_ame_mean)
hist(evapx_ame_mean)
evapx_ame_median<-app(evapx_ame_stack, median)
plot(evapx_ame_median)
hist(evapx_ame_median)
writeRaster(evapx_ame_median, "./WorldClim_30s/wc2.1_30s_vapr/wc2.1_30s_vapr_Ame/wc2.1_30s_vapr_Ame_median.tif", overwrite=TRUE, names="vapr_Ame_median")


#####################################    Radiación solar     ###
#mes mas alto y mes mas bajo
for (x in 1:12) {
  radpx<-rast(paste("./WorldClim_30s/wc2.1_30s_srad/wc2.1_30s_srad_",x,".tif", sep=""))
  radpx_ame<-crop(x=radpx, y=dem_ame)
  writeRaster(radpx_ame, paste("./WorldClim_30s/wc2.1_30s_srad/wc2.1_30s_srad_Ame/wc2.1_30s_srad_",x,"_Ame.tif", sep=""), overwrite=TRUE)
}

rast("./WorldClim_30s/wc2.1_30s_srad/wc2.1_30s_srad_10_MX.tif")

###como estas variables estan mes a mes, entonces debo unir todas las bbdd en un stack y de esta manera entonces calcular un nuevo raster con la media y otro con la mediana para luego escoger
radpx_ame_stack<- rast(list.files(path = "./WorldClim_30s/wc2.1_30s_srad/wc2.1_30s_srad_Ame/", full.names = TRUE))
summary(radpx_ame_stack)

radpx_ame_mean<-app(radpx_ame_stack, mean)
plot(radpx_ame_mean, main="mean")
radpx_ame_median<-app(radpx_ame_stack, median)
plot(radpx_ame_median, main="median")
writeRaster(radpx_ame_mean, "./WorldClim_30s/wc2.1_30s_srad/wc2.1_30s_srad_Ame/wc2.1_30s_srad_Ame_median.tif", overwrite=TRUE, names="srad_MX_median")


# Correlacion -------------------------------------------------------------

###lo que sigue es crear el script para extraer los valores para cada punto y hacer un correlograma para cada especie

### Normalidad?
#Para esto, debo crear un bucle, de tal manera que  el resulatdo final es una matriz donde las filas son las especies y las columnas son las variables y el valor de cada celda es el valor del P de cada variable para cada especie

#cargar las especies
spp.dir<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD", pattern = ".shp$", full.names = T)
#Esto contiene la ruta de cada shapefile


#stack que contiene todas las variables
vbles_stack<-rast(c(list.files(path = "./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/", full.names = TRUE)[-c(18,19,10,11)]))


#aqui va el bucle para la matriz de normalidades
t.sp.norm<- data.frame(
  matrix(nrow = length(spp.dir), ncol = length(vbles_stack@ptr$names)))
names(t.sp.norm)<-vbles_stack@ptr$names
row.names(t.sp.norm)<- list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD", pattern = ".shp$", full.names = F)


#sp.ext está mas abajo sp.ext<-extract(x=vbles_stack, y=points.sp, ID=FALSE)
#necesito estandarizar la variables porque las voy a comparar
for (x in 1:length(spp.dir)) { # 1:length(spp.dir)
  points.sp<-vect(spp.dir[x]) #cargo la especie
  sp.ext.nor<-extract(x=vbles_stack, y=points.sp, ID=FALSE) #extraigo los valores
  sp.ext.nor.sc<-as.data.frame(scale(sp.ext.nor)) #normalizo
  #aplico la prueba
  #problema: Shapiro es entre vectores entre 3 y 5000 datos, debo hacer un if-else
  
  if (length(sp.ext.nor.sc[,1])<5000)
  { shap.list<-apply(sp.ext.nor.sc, MARGIN=2, FUN=shapiro.test)} #se almacena en una lista, donde cada una es el resultado de la prueba
  else {
    sp.ext.nor.sc.spl<-sp.ext.nor.sc[sample(nrow(sp.ext.nor.sc), size=5000), ] #subset de 5000 de tamaño
    shap.list<-apply(sp.ext.nor.sc.spl, MARGIN=2, FUN=shapiro.test)
  }
  for (y in 1:length(shap.list)) {
  t.sp.norm[x,y]<-shap.list[[y]]$p.value
  }
}

length(which(t.sp.norm<0.05)) #1026 son no normales 94.74%
length(which(t.sp.norm>0.05)) #57 son normales 5.26%

#CONCLUSIÓN: correlacion - Asumiré no parametricidad (Spearman)


corrplot.mixed(cor_sp, lower = "number", 
               lower.col ="black", number.cex=.7,
               upper="square", order="FPC", title=" ",
               tl.cex=.7, tl.pos="lt", diag = "u")

cor.matrix<-as.data.frame(cor_sp) #esta función esta modificada en la linea169
#write cada tabla

#con esta función solo se obtienen los valores únicos, osea del triángulo superior
#ahora a armar el bucle para la correlación de las especies

library(corrplot)
#esta funcion extrae los valores de la matriz de correlación y la reordena en 3 columnas
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

spp.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD", pattern = ".shp$", full.names = F) #solo tiene la lista de las especies

#paquetes importantes para el filtro: terra,  corplot, y la función de flttencorrmatrix
for (x in 1:length(spp.dir)) {
  #Abrir shapes
  points.sp<-vect(spp.dir[x]) 
  #extraer valores de los rasters
  sp.ext<-extract(x=vbles_stack, y=points.sp, ID=FALSE)
  #estandarizar valores extraidos
  sp.ext.sc<-as.data.frame(scale(sp.ext))
  #correlación
  cor_sp<- cor(sp.ext.sc, method = "spearman") #sp.ext.sc contiene los valores por pixel de cada punto de presencia estandarizados
      #write matriz
  write.table(cor_sp, paste("C:/Users/Mauricio Diaz/Documents/Maestria/spp. records_PC/specialists_PC/spec.corr/", substr(spp.list[x],1,nchar(spp.list[x])-4), "_corr.txt", sep=""), sep = "\t", dec=".")
  #filtro
  
  cor_sp.long<-flattenCorrMatrix(cor_sp)
  cor_sp.long.filter<-cor_sp.long[which(cor_sp.long[,3]>=-0.8 & cor_sp.long[,3]<=0.8),]
  cor_sp.long.filter$cor<-round(cor_sp.long.filter$cor,3) #con esta función escojo solo los 3 primeros decimales
    #write filtro
# write.table(cor_sp.long.filter, paste("C:/Users/Mauricio Diaz/Documents/Maestria/spp. records_PC/specialists_PC/spec.corr/", substr(spp.list[x],1,nchar(spp.list[x])-4), "_filter.txt", sep=""), sep = "\t", dec=".", row.names=F) 
  #substr(spp.dir[x],1,nchar(spp.dir[x])-4) elimina los últimos 4 caracteres de cada especie (es decir= ".shp")
  
}

####Definición de las M por especie y cortes de variables para cada una
#por ahora voy a hacerlo para una sola especie, luego la idea es hacer un blucle para todas

#de manera manual obtuve la M para A. woodhouseii
A.wod.M<-vect("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.shapes_DD/A. woodhouseii/A. woodhouseii_M.shp")
plot(A.wod.M)

#cortar raster a la M
A.wod.M.stack<-crop(vbles_stack, A.wod.M)
A.wod.M.stack2<-mask(A.wod.M.stack, A.wod.M)
plot(A.wod.M.stack2)

vs<-c(1,7,8,9,10,11,12,15,13,16,19)

#guardar rasters en formato ascii
for (x in 1:length(vs)) {
  writeRaster(A.wod.M.stack2[[x]], 
              paste("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.shapes_DD/A. woodhouseii/M_variables/", names(A.wod.M.stack2)[vs[x]], "_A.wood.asc", sep = ""), NAflag=-9999, overwrite=T)
}


# Corte de variables ------------------------------------------------------

#esta sección tiene el bucle que ava a cortar las variables para cada especie
library(stringr)
spp.dir #direcciones de los registros cada especie
spp.list #nombres del archivo .shp por especie
spp.carp.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = T)
spp.Mshp.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.shapes_DD/", full.names = T, pattern=".shp$")
spp.sets<-read.table("D:/Maestria_DD/Shapes MSc_DD/WorldClim_30s/spec.sets.txt", header=T, sep="\t")


for (x in 1:length(spp.list)) {
  print(spp.list[x])
  #1. Crear carpeta M en la carpeta de cada especie
  dir.create(paste(spp.carp.list[x], "/M_variables", sep = ""))
  #2. Abrir shp de cada especie
  shp.spp <- vect(spp.Mshp.list[x])
  #3. Indentificar variables por especie por set
  subs.spp <- spp.sets[spp.sets[, 1] %in%
                         gsub("_", " ", substr(spp.list[x], 1, (nchar(spp.list[x]) - 4))), ]
  
  #4 Cortar variables
  sets <- levels(factor(subs.spp[, 2])) #cuantos sets hay para la spp x
  for (i in 1:length(sets)) {
    if (length(which(subs.spp[, 2] == sets[i])) >= 1) {
      dir.create(paste(spp.carp.list[x], "/M_variables/", sets[i], sep = "")) #crea vble por set
      #4.1 crear un subset stack donde tenga los raster correspondientes a las variables selec. para el set 1
      subs.stack <-
        vbles_stack[[names(vbles_stack) %in% subs.spp[which(subs.spp[, 2] ==
                                                              sets[i]), 4]]]
      #4.2 cortar el stak a la extensión de la M (shape)
      subs.stack.crop <- crop(x = subs.stack, y = shp.spp)
      subs.stack.mask <- mask(x = subs.stack.crop, mask = shp.spp)
      #4.3 guardar
      writeRaster(subs.stack.mask,
                  paste0(
                    spp.carp.list[x],
                    "/M_variables/",sets[i],"/",names(subs.stack.mask),".asc"),overwrite=T, NAflag= -9999)
    print(sets[i])}
  }
print(spp.list[x])}

# Preparación archivos de especies ----------------------------------------

#Tres archivos csv
  #sp_train: coniene los puntos de entrenamiento para los modelos (80% del total)
  #sp_test: contiene los puntos de testeo para los modelos candidatos (20% del total)
  #sp_joint: contiene los valores completos
spp.carp.list #directorio de las carperas de las especies
spp.list #lista de especies.shp
vbles_stack #stack con las variables
spp.dir #directorio de los shp de las especies

#library(sf)
#library(stringr)
for (x in 1:length(spp.list)) {
  spp.reg<-st_read(spp.dir[x])
  spp.reg.df2<- spp.reg %>%
    st_coordinates() %>%
    as.data.frame()
  spp.reg.df<- spp.reg.df2[,c(1:2)]
  spp.reg.df$species<- paste(
    substr(spp.list[x], 1,1), #obtengo la primera letra
    substr(word(gsub("_", " ", spp.list[x]),start=2,end=2,sep=fixed(" ")),1,4), sep="_")
  names(spp.reg.df)<- c("longitude", "latitude", "species")
  data_joint<- spp.reg.df[,c(3,1,2)]
  
  samples<- sample(2, nrow(data_joint), replace=T, prob = c(0.8,0.2))
  data_train<- data_joint[samples == 1,]
  data_test<- data_joint[samples == 2,]
    
  write.csv(data_train, paste0(spp.carp.list[x], "/", "Sp_train.csv"), row.names = F)
  write.csv(data_test, paste0(spp.carp.list[x], "/", "Sp_test.csv"), row.names = F)
  write.csv(data_joint, paste0(spp.carp.list[x], "/", "Sp_joint.csv"), row.names = F)
}

#Fin del script
     