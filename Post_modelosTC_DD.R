library(biscale)
library(ggpmisc)
library(ggridges)
library(janitor)
library(patchwork)
library(sf)
library(stringr)
library(terra)
library(tidyverse)

#En este script voy a hacer todo lo que lleva el postprocesamiento de los modelos   

# ------------------------------------- Binarización de los modelos --------------------------------


mods.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = T) #modelos direcciones 
spp.mod.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F)

tabla.var.bin<-data.frame(spp=list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F), `10ptp.median`=NA, `10ptp.mean`=NA, `10ptp.max`=NA, `10ptp.rep.max`=NA, `10ptp.min`=NA, `10ptp.rep.min`
                          =NA, Mean_AUC_ratio=NA, pval_pROC=NA, Omiss.rat=NA, AICc=NA, delta_AICc=NA, W_AICc=NA, num_par=NA)

#max.res[,39] #X10.percentile.training.presence.Logistic.threshold 

#Primera parte: Extracción de 10 percentil training presence por especie y algunas métricas
#comparador
tab.mods<-read.table("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.lista.txt", sep="\t", dec=".", header=T)


for (x in 1:length(tab.mods[,1])){
  #1. arbir carpeta con modelos finales
  carp.fm<- list.files(paste0(mods.dir[x], "/Final_Models"), full.names=F)
  #1.1 Aplicar if else porque el órden de las carpetas es diferente al nombre de las variables
  if (tab.mods[x,1]==substr(spp.mod.list,1,nchar(spp.mod.list)-2)[x]) { 
    #2. Calculo de la mediana
    #2.1 abrir maxent results
    max.res<-read.csv(paste0(paste0(mods.dir[x], "/Final_Models/", tab.mods[x,2], "/maxentResults.csv")))
    #2.2 Calculo de estadisticos para el 10PTP LT
    tp.tp.lh<-max.res$X10.percentile.training.presence.Logistic.threshold
    tabla.var.bin[x,2]<-median(tp.tp.lh[1:10]) #median
    tabla.var.bin[x,3]<-mean(tp.tp.lh[1:10]) #median
    tabla.var.bin[x,4]<-max(tp.tp.lh[1:10]) #max
    tabla.var.bin[x,5]<-max.res[which(
      max.res$X10.percentile.training.presence.Logistic.threshold==max(tp.tp.lh[1:10])),1] #max.rep
    tabla.var.bin[x,6]<-min(tp.tp.lh[1:10]) #min
    tabla.var.bin[x,7]<-max.res[which(
      max.res$X10.percentile.training.presence.Logistic.threshold==min(tp.tp.lh[1:10])),1] #min.rep
    #2.3 calculos de resultados de calibración
    cal.res<-read.csv(paste0(paste0(mods.dir[x], "/Calibration_results/calibration_results.csv")))
    tabla.var.bin[x,8]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,2]
    tabla.var.bin[x,9]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,3]
    tabla.var.bin[x,10]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,4]
    tabla.var.bin[x,11]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,5]
    tabla.var.bin[x,12]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,6]
    tabla.var.bin[x,13]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,7]
    tabla.var.bin[x,14]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,8]
  }
}

#write.csv2(tabla.var.bin, "F:/Maestria_DD/spp.records_DD/Info.FM.csv", row.names = F, dec=".")
tabla.var.bin<-read.csv2("F:/Maestria_DD/spp.records_DD/Info.FM.csv")

#Segunda parte: Binarizar los modelos usando el umbral de cada especie


for (x in 1:8) { #length(tab.mods[,1])
  #1. Abrir el raster por especie
  ras.median<-rast(paste0(mods.dir[x], "/Final_Models/", tab.mods[x,2], "/",
                          paste(substr(spp.mod.list, 1,1), #obtengo la primera letra
                                substr(word(gsub("_", " ", spp.mod.list),start=2,end=2,sep=fixed(" ")),1,4), sep="_")[x], 
                          "_median.asc"))
  ras.median2<-ras.median
  ras.median2[ras.median2>tabla.var.bin[x,2]]<-1
  ras.median2[ras.median2<=tabla.var.bin[x,2]]<-0
  writeRaster(ras.median2, paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_raw/", 
                                  spp.mod.list[x] |> #nombre del .tif de cada especie
                                    word(1, sep = fixed("_")) |> #separo la primera palabra (genero)
                                    str_sub(1,2) |> #solo extraigo la primera letra
                                    paste( #creo el paste para unir genero y especie
                                      spp.mod.list[x] |> #nombre del .tif de cada especie
                                        word(2, sep = fixed("_")) |> #extraigo la segunda palabra
                                        substr(1,4), sep="_")#extraigo los 4 primeros caracteres
                                  , "_bm.tif"), overwrite=TRUE) #binarian model
  print(spp.mod.list[x])}
#-------------------------------------------------------------------------------------------------------------- # 

# -------------------------------------------------Mapas de riqueza ------------------------------------------

#-------------------------------------------------------------------------------------------------------------- # 

# Raster de riqueza total -------------------------------------------------

spp.bin.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", pattern=".tif$", full.names = T) #ubicación de la carpeta con los Modelos_binarios
spp.bin.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", pattern=".tif$", full.names = F)
ext.ref<-rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/wc2.1_30s_bio_1_Ame.tif") #raster con la extension de referencia con la que quiero cortar los raster de distribución 

# abro todos los rasters en una lista
bin.list<-list()
for (i in 1:length(spp.bin.dir)) {
  bin.list[[i]]<-rast(spp.bin.dir[i])
  names(bin.list)[i]<-substr(spp.bin.list[i], 1,nchar(spp.bin.list[i])-4)
}

#calculo de riqueza
bin.list |> 
  lapply(FUN=extend, y=ext(ext.ref)) |> 
  rast() |> 
  sum(na.rm = T) |> 
  writeRaster("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_total.tif", 
              overwrite=T)

#1. Recortar distribuciones a México ------------------------------------------------- 

#raster con la extension de referencia con la que quiero cortar los raster de distribución 
ext.ref.mx<-vect("F:/Maestria_DD/Shapes_MSc_DD/Mexico_Estados/Mexico_continent.shp")
spp.bin.dir #direccion de los Modelos_binarios
spp.bin.list #lista de nombres de modelos binarios

for (x in 1:length(spp.bin.dir)) {
  spp.bin.dir[x] |>
    rast() |>
    mask(ext.ref.mx, touches=F) |>
    crop(ext.ref.mx) |> 
    writeRaster(
      paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/",
             substr(spp.bin.list[x], 1, nchar(spp.bin.list[x]) - 4),
             "_MX.tif"), overwrite=T)
  print(spp.bin.list[x])
}

#2. Calculo de riqueza en México -------------------------------------------------------

bin.list.mx<-list()
spp.bin.dir.mx<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = T) # Ubicación de la carpeta con los Modelos_binarios
spp.bin.list.mx<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = F) 

# Abro todos los rasters en una lista
for (i in 1:length(spp.bin.dir.mx)) {
  bin.list.mx[[i]]<-rast(spp.bin.dir.mx[i])
  names(bin.list.mx)[i]<-substr(spp.bin.list.mx[i], 1,nchar(spp.bin.list.mx[i])-4)
}

# Cálculo de riqueza
bin.list.mx |> 
  lapply(FUN=extend, y=ext(ext.ref.mx)) |> 
  rast() |> 
  sum(na.rm=T) |> 
  writeRaster("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif", overwrite=T) #|> plot()


# ----------------------------------------------- T y C S7 ----------------------------------------------------  
# Corte a categorias ------------------------------------------------------

# Direcciones de los shapes de cambios
shp.cambios.dir<- list.files("F:/Maestria_DD/Shapes_MSc_DD/INEGI/Cambios/TC_s7/", pattern = ".shp$", full.names = T)
# Lista de los cambios "T - C"
shp.cambios.list<-list.files("F:/Maestria_DD/Shapes_MSc_DD/INEGI/Cambios/TC_s7/", pattern = ".shp$", full.names = F)
shp.cambios.list<-substr(shp.cambios.list, 1,1) #eliminar _s7

# Ubicación de direcciones de las especies cortadas a México
spp.bin.dir.mx<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = T) 
# Lista de los nombres de los archivos de las especies
spp.bin.list.mx<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = F)

# abro todos los rasters en una lista
bin.list.mx<-list()
for (i in 1:length(spp.bin.dir.mx)) {
  bin.list.mx[[i]]<-rast(spp.bin.dir.mx[i])
  names(bin.list.mx)[i]<-substr(spp.bin.list.mx[i], 1,nchar(spp.bin.list.mx[i])-4)
}


for (x in 1:length(shp.cambios.list)) {# contador de las categorias
  print(paste0("-----", shp.cambios.list[x], "-----"))
  cat.shp<-vect(shp.cambios.dir[x])
  for (y in 1:length(bin.list.mx)) {
    bin.list.mx[[y]] |> 
      mask(mask=cat.shp, touches=F) |> 
      writeRaster(paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/Z_",
                         substr(shp.cambios.list[x],1,2),"/", spp.bin.list.mx[y] |> #nombre del .tif de cada especie
                           word(1, sep = fixed("_")) |> #separo la primera palabra (genero)
                           str_sub(1,2) |> #solo extraigo las dos primeras letras
                           paste( #creo el paste para unir genero y especie
                             spp.bin.list.mx[y] |> #nombre del .tif de cada especie
                               word(2, sep = fixed("_")) |> #extraigo la segunda palabra
                               substr(1,4),sep="_"), #extraigo los 4 primeros caracteres
                         "_",substr(shp.cambios.list[x],1,2),".tif"), overwrite=T)
    print(spp.bin.list.mx[y])}
}


# Calculo de riqueza por cat ----------------------------------------------
cat.dir.TC<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", full.names = T, pattern=".tif$") # dirección de las carpetas de categorias
cat.list.TC<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/",full.names = F,pattern=".tif$")#lista de categorias
ext.ref.mx #extend referencia de mexico

# abro todos los rasters en una lista
spp.list.TC<-list()
for (x in 1:length(cat.dir.TC)) { #x=contador de categorias
  for(y in 1:56) { #y=contador de las especies
    #1. Abrir raster en una lista
    spp.list.TC[[y]]<-rast(list.files(cat.dir.TC[x], full.names = T, pattern = ".tif$")[y])
  }
  spp.list.TC |>
    lapply(FUN = extend, y = ext(ext.ref.mx)) |>
    rast() |>
    sum(na.rm = T) |>
    writeRaster(
      paste0(
        "F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/",
        cat.list.TC[x],
        ".tif"),
      overwrite = T)
}


# Pixeles de cada especie por categoría ---------------------
cat.dir.TC # dirección de las carpetas de categorias
cat.list.TC #lista de categorias
ext.ref.mx #extend referencia de mexico

spp.cat.area<-data.frame(spp=list.files(cat.dir.TC[1]) |> #lista de especies
                           str_sub(1, -8),
                         C=NA, #categorias
                         T=NA,
                         Tot.dist=NA) #Pixeles totales de la distribución

list.spp.pix<-list()
for (x in 1:length(cat.dir.TC)) {
  print(paste0("-----",str_sub(cat.list.TC[x],3,4),"-----"))
  list.spp<-list() #lista que contendra los DF por categorias
  for (y in 1:length(list.files(cat.dir.TC[x]))) { 
    print(str_sub(list.files(cat.dir.TC[x])[y], 1, -7))
    #1. Abrir los raster
    list.spp[[y]]<-rast(list.files(cat.dir.TC[x], full.names = T, pattern = ".tif$")[y])
    names(list.spp)[y]<-str_sub(list.files(cat.dir.TC[x])[y], 1, -7)
  }
  list.spp.pix[[x]]<- list.spp |> 
    lapply(FUN=extend, y=ext(ext.ref.mx)) |> #Modificar el extend de los rasters
    rast() |> 
    freq() |> #calcular la frecuencia
    mutate(cat=str_sub(cat.list.TC[x],3,4)) |> #Crear una columna con el nombre de cada categoría
    rename(spp=1, #cambio nombres de columnas
           value=2,
           n.pix=3,
           cat=4)
  names(list.spp.pix)[x]<- str_sub(cat.list.TC[x],3,4)
}

#converitr lista en una tabla
list.spp.pix.P<- as.data.frame(do.call(rbind, list.spp.pix))

#cambiar los valores por los nombres de las especies
for (y in 1:56) {
  list.spp.pix.P[which(list.spp.pix.P[,1]==y),1] <- str_sub(list.files(cat.dir.TC[1]),1,-7)[y]
}


# --------------------------------- Familia y ordenes  ---------------------------------------
spp.nom.cods<-read.table("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/spp.famord.txt", sep="\t", dec=".", header = T)
NACC<-read.csv("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/NACC_list_species.csv")

order<-NACC |> 
  select(order) |> 
  distinct() |> 
  mutate(count=1:31)

family<-NACC |> 
  select(family) |> 
  distinct() |> 
  mutate(count=1:133)

spp.nom.cods$orden.cod<-NA
spp.nom.cods$famil.cod<-NA
names(spp.nom.cods)

for (x in 1:nrow(spp.nom.cods)) {
  
  for (y in 1:nrow(family)) {
    if(spp.nom.cods[x,3]==family[y,1]) {spp.nom.cods[x,7]<-family[y,2]}
  }
  
  for (z in 1:nrow(order)) {
    if(spp.nom.cods[x,4]==order[z,1]) {spp.nom.cods[x,6]<-order[z,2]}
  }
}

############################################## Columnas numero total de pixeles
#1. en méxico
spp.bin.dir.mx #direccion de las especies cortadas a Mexico
spp.bin.list.mx #lista de las especies

#distribuciones para mexico
list.dist<-list()
for (x in 1:length(spp.bin.dir.mx)) {
  #1. Abrir rasters en una lista
  list.dist[[x]] <- rast(spp.bin.dir.mx[x])
  names(list.dist)[x] <- str_sub(spp.bin.list.mx[x], 1,-11)
}


#guardar tabla con las distribuciones
tab.distxcats<- list.spp.pix.P |>
  #1. unir la tabla con los nombres a traves del codigo de la especie
  merge(y = spp.nom.cods,
        by = "spp") |> 
  #2. seleccionar solo los pixeles de presencia
  filter(value == 1) |> 
  #cambiar formato a ancho
  pivot_wider(names_from = cat,
              values_from = n.pix) |>
  #2.1 cambiar posicion de la columna nom
  relocate(name, .after = spp) |>
  #3. pixeles de distribucion total
  mutate(n.pix.MX = list.dist |>
           #3. 1 Conteo total de pixeles de distribucion en México para cada spp
           lapply(FUN = extend, y = ext(ext.ref.mx)) |> 
           rast() |>
           freq() |> 
           #3.1 seleccionar solo los pixeles de presencia
           filter(value==1) |> 
           select(count) |> 
           pull(),
         #4. porcentaje de pixeles de distribucion de cada especie en cada categoria
         pc.dist.C= C*100/list.spp.pix.P |> 
           filter(value==1) |> 
           group_by(spp) |> 
           summarise(count=sum(n.pix)) |> 
           select(count) |> 
           pull(), 
         pc.dist.T= T*100/list.spp.pix.P |> 
           filter(value==1) |> 
           group_by(spp) |> 
           summarise(count=sum(n.pix)) |> 
           select(count) |> 
           pull()) 

# tab.distxcats |> write.table("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/tab.distxcats.txt", sep="\t", dec=".", row.names=F)

tab.distxcats<- read.table("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/tab.distxcats.txt", sep="\t", dec=".", header=T)
head(tab.distxcats)

sum(tab.distxcats[1,12:13]) #100% ok

tab.distxcats[which(tab.distxcats$endemismo=="E"),3]<-1
tab.distxcats[which(tab.distxcats$endemismo=="Q"),3]<-2
tab.distxcats[which(tab.distxcats$endemismo=="NE"),3]<-3

tab.distxcats |> 
  pivot_longer(cols=c(pc.dist.C, pc.dist.T),
               names_to= "pc.cats",
               values_to = "pc") |> 
   # group_by(nom, pc.cats) |> 
   # summarise(count=sum(pc)) |> 
  ggplot(aes(x=fct_reorder(name, Posicion,.desc=F), y=pc, fill=pc.cats))+
  geom_bar(stat="identity") +
  labs(y="Proportion", fill="Categories") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 13),
        axis.text.y = element_text(size = 13),
        panel.background = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = 50, col="black", linetype="dashed")
  #scale_fill_manual(values=c("#248f5d","#e5b636"), labels = c("Conserved","Transformed"))

ggsave(file = "./outputs/catsxspp.svg",
       width = 1706,
       height = 1271,
       scale=3,
       units ="px")

tab.distxcats |> 
  pivot_longer(cols=c(pc.dist.C, pc.dist.T),
               names_to= "pc.cats",
               values_to = "pc") |> 
  filter(pc>75) |> 
  select(name) |> 
  pull()


#proporción de las especies?
tab.distxcats |> 
  filter(spp=="Tr_ambi") #spp=="Gl_hosk" | spp=="Ca_oliv" | spp=="Ca_mela"
mean(tab.distxcats$pc.dist.C)
mean(tab.distxcats$pc.dist.T)

#cual especie tiene mas distribucion en T que en C?
tab.distxcats |> 
  mutate(rest=pc.dist.C-pc.dist.T) |> 
  arrange(desc(rest)) |> 
  head()

50*100/56 #especies con la mayoría de su distribucion en T
6*100/56 #especies con la mayoría de su distribucion en C

#---------------------------------------- Relación con la pendiente -----------------------------------------------
mx.slope<- rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
mx.slope |> 
  as.data.frame() |> 
#  na.omit() |> 
  nrow()

TC.Z.dir<- list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", pattern=".tif$", full.names = T) # Direcciones de la riqueza por categoría
TC.Z.list<- list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", pattern=".tif$", full.names = F) # Lista de los nombres de los rasters de riqueza por categoría
cat.dir.TC


bl.stack.cat.TC<- rast(TC.Z.dir)
#1. convertir raster de categorías a puntos

cat.Z.list.TC<-list()
for(i in 1:length(cat.dir.TC)){
  # 1. convertir raster de riqueza por categorías en puntos
  cat.point<-as.points(bl.stack.cat.TC[[i]]) #contiene los rasters de riqueza por categoria
  # 2 Convertir tabla en un DF y guardarlo en una lista
  cat.Z.list.TC[[i]]<-
    # 2.1 Extraer los valores de pendiente por valor de riqueza
    terra::extract(mx.slope, cat.point, bind=T, xy=T) |> 
    as.data.frame() |> 
    # 2.2 Agregar columna de categoría
    mutate(cat= cat.list.TC[i] |> 
             str_sub(3,3)) |> 
    # 2.3 Renombrar variables
   dplyr::rename(z=sum, slope=slope_mx_g) |> 
    # 2.4 Reordenar
    relocate(c(x,y)) |>  
    # 2.5 filtrar
    filter(z!=0)
    # 2.6 cambiar los nombres
  names(cat.Z.list.TC)[i]<-paste0(substr(cat.list.TC[i], 3,3),"_Z")
  # 3 guardar el txt
  write.table(cat.Z.list.TC[[i]], paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", levels(factor(cat.Z.list.TC[[i]]$cat)), "_Z_slp.txt"), sep="\t", dec=".", row.names = F)
}

#unirlos todos en un solo DF
Z.slp.cats<-as.data.frame(do.call(rbind, cat.Z.list.TC)) #contiene todos los valores de pendiente por valor de riqueza en cada categoria
#Z.slp.cats$count<-1 #agrego una columan de contador

#Z.slp.cats<-
  Z.slp.cats |> 
  #elimino NA
  na.omit() #|> write.table("./outputs/Z.slp.cats.txt", sep = "\t", dec=".", row.names = F)
summary(Z.slp.cats) #verifico

# valores de pendiente por pixel de riqueza
Z.slp.cats<- read.table("F:/Maestria_DD/Shapes_MSc_DD/outputs/tablas/Z.slp.cats.txt", sep = "\t", dec=".", header=T)
head(Z.slp.cats)

#Cuál es la mediana de la pendiente por categoria?

#options(pillar.sigfig =4)
Z.slp.cats |> 
  group_by(cat) |> 
  # summarize(median=median(slope),
  #           sd=sd(slope))
  summarise(cat.z=n()) |> 
  mutate(tot.pix=(398888100-396407315),
         cat.pix=c(1236896, 1243889),
         cat.sinz=cat.pix-cat.z,
         #cuántos pixeles de T y C estan ocupados por al menos una especie en proporción
         p.z=cat.z*100/cat.pix, 
         p.sinz.tot=cat.sinz*100/cat.pix) |> 
  relocate(c(tot.pix, cat.pix), .before=cat.z) |> 
  pivot_longer(c(p.z, p.sinz.tot), names_to="props", values_to = "values") |> 
#  mutate(props2=paste0(props,"_",cat)) |> 
  ggplot(aes(x=cat,y=values, fill=props))+
    geom_bar(stat="identity") +
  labs(y="Percentage")+
  scale_x_discrete(labels=c("Conserved","Transformed"))+
#  scale_fill_discrete(name = " ", labels = c("Without species", "With species"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")


ggsave(filename = "./outputs/richness_cat2.svg",
         width = 10,
         height = 10, #alto
         scale=1,
         units ="cm",
         dpi = 200)


# densidades de kernel ----------------------------------------------------

# Slope -------------------------------------------------------------------
#plot total
slp_dens_tot<- 
  ggplot(Z.slp.cats, aes(x=slope, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#e5b636"),labels=c("Conserved", "Transformed")) + 
  scale_color_manual(values=c("#248f5d","#e5b636"), labels=c("Conserved", "Transformed")) +
  labs(x="Terrain slope", y="Density", fill="Categories", color="Categories") +
  theme_classic() +
  theme(legend.position="none")


# ggsave(filename = "./outputs/slp_dens_tot.png",
#          width = 15,
#          height = 10, #alto
#          scale=2,
#          units ="cm",
#          dpi = 200)



#por cada uno
f_labels <- data.frame(cat = c("C", "T"), label = c("A", "B"))
#f_slop<- #esto es para el plot doble
ggplot(Z.slp.cats, aes(x=slope, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#e5b636")) + 
  scale_color_manual(values=c("#248f5d","#e5b636")) +
  facet_wrap("cat", labeller = labeller(cat=c(C="Conserved",
                                              T="Transformed"))) +
  theme_classic()+
  theme(legend.position="none",
        strip.text =  element_text(size=10, face="bold"))+
  labs(x="Terrain slope", y="Density")
  #agregar labels (A - B)
 #+ geom_text(x = 33, y = 0.53, aes(label = label), data = f_labels, color="black", fontface = "bold") #hay que crear un objeto para esto: f_labels

# ggsave(filename = "./outputs/riq_dens_facet.png",
#          width = 15,
#          height = 10, #alto
#          scale=2,
#          units ="cm",
#          dpi = 200)


# Densidad riqueza -----------------------------------------------------------------
#plot total
ggplot(Z.slp.cats, aes(x=z, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#e5b636"),labels=c("Conserved", "Transformed")) + 
  scale_color_manual(values=c("#248f5d","#e5b636"), labels=c("Conserved", "Transformed")) +
  labs(x="Richness", y="Density", fill="Categories", color="Categories")

# ggsave(filename = "./outputs/riq_dens_tot.png",
#          width = 15,
#          height = 10, #alto
#          scale=2,
#          units ="cm",
#          dpi = 200)

#Densidad por cada uno
f_labels
f_Z<- #esto es para el plot doble
ggplot(Z.slp.cats, aes(x=z, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#e5b636")) + 
  scale_color_manual(values=c("#248f5d","#e5b636")) +
  # facet_wrap("cat", labeller = labeller(cat=c(C= "Conserved",
  #                                             T= "Transformed"))) +
  theme_classic()+
  labs(x="Species richness", y="Density")+
  theme(legend.position="none",
        strip.text =  element_text(size=10, face="bold")) 
#agregar labels (A - B)
#+ geom_text(x = 29, y = 0.26, aes(label = label), data = f_labels, color="black", fontface = "bold") #hay que crear un objeto para esto: f_labels

# ggsave(filename = "./outputs/riq_dens_facet.png",
#          width = 15,
#          height = 10, #alto
#          scale=2,
#          units ="cm",
#          dpi = 200)

# library(patchwork)
# f_Z/f_slop +
#   plot_annotation(tag_levels = "A")

# ggsave(filename = "./outputs/riqslope_dens_facet.png",
#        width = 10,
#        height = 15, #alto
#        scale=2,
#        units ="cm",
#        dpi = 200)

# Plots de Riqueza vs pendiente --------------------------------------------------------



#grafico total 
# densidades
Z.slp.cats |> 
  sample_n(size=100) |> 
  ggplot(aes(y=factor(z), x=slope, width = after_stat(density))) + 
  geom_density_ridges(scale = 0.1)+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5)+ #plotee la mediana dentro de cada uno
  geom_vline(xintercept= 2.318058, #esta es la mediana de los pixeles de la riqueza de presencia potencial
              color="red",linetype = "dashed") +
  geom_vline(xintercept= c(0,2,10,20),
             color="black") +
  labs(y="Species richness", x="Terrain slope") +
  facet_wrap(~cat)+
  theme_classic()

ggsave(filename = "./outputs/slp_bxplt2.svg",
       width = 15,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)

# boxplot
Z.slp.cats |> 
#  sample_n(size=100) |> 
  ggplot(aes(x=slope, y=factor(z))) + 
  geom_boxplot() +
  geom_vline(xintercept= 2.318058, #esta es la mediana de los pixeles de la riqueza de presencia potencial
             color="red",linetype = "dashed")+
  geom_vline(xintercept= c(0,2,10,20),
             color="black") +
  labs(y="Species richness", x="Terrain slope") +
  scale_x_continuous(trans="pseudo_log", limits=c(0, 35))+
  theme_classic()

ggsave(filename = "./outputs/slp_bxplt.svg",
       width = 8,
       height = 10, #alto
       scale=1.5,
       units ="cm",
       dpi = 200)

Z.slp.cats |> 
#  sample_n(size=100) |> 
  filter(cat=="C") |> 
  ggplot(aes(x=slope, y=factor(z), fill=cat)) + 
  geom_boxplot() +
  geom_vline(xintercept= Z.slp.cats |> 
               filter(cat=="C") |> 
               summarise(median=median(slope)) |> 
               pull(), #esta es la mediana de los pixeles de la riqueza de presencia potencial en zonas conservadas 3.38
             color="red",linetype = "dashed")+
  geom_vline(xintercept= c(0,2,10,20),
             color="black") +
  labs(y="Species richness", x="Terrain slope") +
#  facet_wrap(~cat, labeller=as_labeller(c(C="Conserved areas",T="Transformed areas")))+
  scale_fill_manual(values=c("#248f5d"))+
  scale_x_continuous(trans="pseudo_log", limits=c(0, 35))+
  theme_classic()+
  theme(legend.position = "none")
Z.slp.cats |> 
  #  sample_n(size=100) |> 
  filter(cat=="C") |> 
  summarise(mean=mean(slope)) |> 
  pull()

ggsave(filename = "./outputs/z_slp_bxplt_C.svg",
       width = 20,
       height = 14, #alto
       scale=1,
       units ="cm",
       dpi = 200)


# mediana de la pendiente por cat -----------------------------------------------------------------
slp_z_med<-
  Z.slp.cats |> 
  group_by(cat, z) |> 
  summarise(median= median(slope)) |> 
  ggplot(aes(y=z, x=median))+
  geom_point()+
  geom_smooth(method="lm", formula=y~x, color="black")+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  labs(x="Terrain slope (median)", y="Species richness") +
  theme_classic()



#slp_z_med2<-
  Z.slp.cats |> 
  group_by(cat, z) |> 
  summarise(median= median(slope)) |> 
  ggplot(aes(y=z, x=median, color=cat))+
  geom_point()+
  #geom_smooth(method="lm", formula=y~x, se=F)+
  geom_abline(aes(intercept=7.02377, slope=1.47603), color="#248f5d") +
  geom_abline(aes(intercept=5.99293, slope=1.76647), color="#e5b636") +
  # stat_poly_eq(formula = y~x, 
  #              aes(label = paste(after_stat(eq.label), ..rr.label.., sep = "~~~")), 
  #              parse = TRUE)+
  scale_color_manual(values=c("#248f5d","#e5b636"), labels = c("Conserved","Transformed"))+
  labs(x="Terrain slope (median)", y="Species richness", color="Categories") +
  theme_classic()+
  theme(axis.title.y = element_blank())

  

slp_z_med+slp_z_med2+ plot_annotation(tag_levels = 'A')

ggsave(filename = "./outputs/slp_x_z.png",
       width = 20,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)



(slp_z_med / slp_spp_bxplt_facet3) + plot_annotation(tag_levels = 'A')

ggsave(filename = "./outputs/slp_merge3.png",
       width = 10,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)

#Las especies se concentran en valores más altos de pendiente


# Slope x categoria -------------------------------------------------------

Z.slp.cats |> 
#  sample_n(size=100) |> 
  ggplot(aes(x=cat, y=slope, fill=cat))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  labs(y="Pendiente del terreno") +
  scale_fill_manual(values=c("#248f5d","#e5b636"))+
  scale_y_continuous(trans="pseudo_log", limits=c(0, 35))+
  geom_hline(yintercept = c(0,2,10,20), color="black") +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y =element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
#mean(Z.slp.cats$slope)

ggsave(filename = "./outputs/slp_cat2.svg",
       width = 10,
       height = 10, #alto
       scale=1,
       units ="cm",
       dpi = 200)

spp_slp.join |> 
  summarise(medi=median(slope, na.rm = T), sd=sd(slope), .by=cat)

# Zonificación de la pendiente en T y C ------------------------------------

reg.bio<-vect("Provincias biogeograficas/Regiones_biog/Regiones_biog.shp")
biomas<-vect("Provincias biogeograficas/Ecoregions2017_dinerstein/Biomas_2017_Dinerstein.shp")
cat.rast<-rast("INEGI/Cambios/INEGI_VII_TC.tif")
z.mx<-rast("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif")
mx.slope<- rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
plot(cat.rast)

#esto es para las regiones o para los biomas dependiendo de que quiera sacar, para no repetir el código que es el mismo
zon.list<-list()
for(i in 1:length(biomas)){
  zon.list[[i]]<-
    terra::extract(cat.rast,
                   terra::extract(
                     mx.slope,
                     z.mx |> crop(biomas[i, ], mask = T) |> as.points(),
                     ID = F,
                     bind = T),
                   ID=F, 
                   bind=T) |> 
    as.data.frame() |> 
    mutate(zonif=biomas[i,]$BIOME_NAME)
  print(biomas[i,]$BIOME_NAME)
}

as.data.frame(do.call(rbind, zon.list)) |> 
  na.omit() |> 
  rename(slope=slope_mx_g, z=sum) |> 
  mutate(cat=recode(INEGI_VII_TC, `1`="C", `2`="T")) |> 
  write.table("./outputs/tablas/zon.biom.txt", sep="\t", dec=".", row.names=F)

zon<-read.table("./outputs/tablas/zon.biom.txt", sep="\t", dec=".", header=T) #ojo debe cambiar el ojbeto dependiendo de la zonificacion
head(zon)

# Kruskal-wallis test -------------------------------------------------------------------
# esto se hace en línea con lo sugerido por Juli, donde me hace notar que es bueno evaluar diferencias significativas dentro de cada zonificiación

zon.lev<-levels(factor(zon$zonif))
aov.table<-data.frame(zon=NA, p.val=NA)
boots<-rep(100,500)
aov.list<-list()

for (x in 1:length(zon.lev)) {

  for(i in 1:length(boots)){
    anova<-kruskal.test(slope ~ cat, data=zon |> 
                          filter(zonif==zon.lev[x]) |> 
                          sample_n(size=boots[i], replace=T))
    
    aov.table[i,1]<-zon.lev[x]
    aov.table[i,2]<-format(anova$p.value, scientific = FALSE) |> 
      as.numeric() |> 
      round(digits=4)
  }
  aov.list[[x]]<-aov.table
}

aov.df<-as.data.frame(do.call(rbind, aov.list))

aov.df |> 
  summarise(p.val.median=median(p.val), p.val.mean=mean(p.val), .by=zon)

aov.df |> 
  ggplot(aes(x=p.val))+
  geom_histogram()+
  geom_vline(xintercept = 0.05)+
  facet_wrap(~zon)

zon |> 
#  sample_n(size=200) |> 
  filter(z!=0) |> 
  filter(zonif!="Mangroves") |> 
  ggplot(aes(x=cat, y=slope, fill=cat))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values=c("#248f5d","#e5b636"))+
  facet_wrap(~zonif, scales="free")+
  scale_y_continuous(trans="pseudo_log", limits=c(0, 35))+
   geom_hline(yintercept = c(0,2,10,20), color="black") +
  theme_classic()+
  labs(y="Pendiente del terreno")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "./outputs/slp_cat_biom.svg",
       width = 22, #reg es 15 x 10
       height = 15, #alto 
       scale=1,
       units ="cm",
       dpi = 200)

# Valores de pendiente dentro de C y T dentro de cada bioma
zon |> 
  #  sample_n(size=200) |> 
  filter(z!=0) |> 
  filter(zonif!="Mangroves") |> 
  summarise(median.slope=median(slope),sd.slope=sd(slope), n=n(), .by=c(zonif, cat)) |> 
  pivot_wider(names_from=cat, values_from = c(median.slope, sd.slope,n))
  ggplot(aes(x=cat, y=median.slope))+
  geom_bar(stat='identity')+
  facet_wrap(~zonif, scales="free")

#construccion de los kruskal para el total
  
  aov.table2<-data.frame(zon=NA, p.val=NA)
  
  for(i in 1:length(boots)){
    anova<-kruskal.test(slope ~ cat, data=zon |> 
                          sample_n(size=boots[i], replace=T))
    
    aov.table2[i,1]<-"total"
    aov.table2[i,2]<-format(anova$p.value, scientific = FALSE) |> 
      as.numeric() |> 
      round(digits=4)
  }
  
  aov.table2|> 
    summarise(p.val.median=median(p.val), p.val.mean=mean(p.val))
  
  aov.table2 |> 
    ggplot(aes(x=p.val))+
    geom_histogram()
  
# Pendiente total ---------------------------------------------------------

zon |> 
    ggplot(aes(x=cat, y=slope, fill=cat))+
    geom_violin()+
    geom_boxplot(width=0.1)+
    scale_fill_manual(values=c("#248f5d","#e5b636"))+
    scale_y_continuous(trans="pseudo_log", limits=c(0, 35))+
    geom_hline(yintercept = c(0,2,10,20), color="black")+
    theme_classic()+
    labs(y="Pendiente del terreno")+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  ggsave(filename = "./outputs/slp_cat3.svg",
         width = 10,
         height = 10, #alto
         scale=1,
         units ="cm",
         dpi = 200)
  
  
# ------------------------------------------- Pendente por especie ----------------------------

cat.dir.TC # dirección de las carpetas de categorias
cat.list.TC #lista de categorias

list.cat.slp<-list()
for(x in 1:length(cat.dir.TC)){
  print(cat.list.TC[x])
  list.spp<-list()
  for(y in 1:56){
    print(list.files(cat.dir.TC[x], full.names=F)[y])
    #1. Abrir el raster de cada especie y convertir a puntos
    cat.point.spp<- rast(list.files(cat.dir.TC[x], full.names=T)[y]) |> 
      as.points()
    #2. Extrater los valores de pendiente por cada pixel de presencia
    tabla<-as.data.frame(terra::extract(x=mx.slope, y=cat.point.spp, bind=T)) #contienes valores de cero
    #3. ordenar y guardar en un objeto de la lista
    list.spp[[y]]<- tabla |> 
      mutate(cat=substr(cat.list.TC[x], 3,4)) |> 
      #3.1 agregar columna con la categoria
      rename(spp=1, 
      #3.2 renombrar columnas
             slope=2) |>
      subset(spp!=0) |> 
      #3.3 solo aquellas con valores diferentes de 0
      mutate(spp=substr(list.files(cat.dir.TC[x])[y],1,nchar(list.files(cat.dir.TC[x])[y])-6)) #3.4 agregue nombre de la spp
  }
  #4. convertir la lista en un Df
  list.cat.slp[[x]]<-as.data.frame(do.call(rbind, list.spp))
  names(list.cat.slp)[x]<- cat.list.TC[x]
  print(cat.list.TC[x])
}

as.data.frame(do.call(rbind, list.cat.slp)) #|>  write.table("./outputs/tablas/spp_slp.txt", dec=".", sep="\t", row.names = F)
# rm(list.cat.slp)

#Barplot
spp_slp<-read.table("./outputs/tablas/spp_slp.txt", dec=".", sep="\t", header=T)
head(spp_slp)
head(spp.nom.cods)

#para ordenar por codigo cada familia y orden debe tener un numero
fams.cods<- data.frame(familia=levels(factor(spp.nom.cods$familia)), cod=1:26)
ord.cods<-data.frame(orden=levels(factor(spp.nom.cods$orden)), cod=1:10)

spp.nom.cods$o.cods<-NA #o.cods= ordenes;  f.cods= familias
spp.nom.cods$f.cods<-NA #o.cods= ordenes;  f.cods= familias

for(x in 1:26){
  for(y in 1:56){
    if(spp.nom.cods$familia[y]==fams.cods$familia[x]) 
      spp.nom.cods$f.cods[y] <- fams.cods$cod[x]
  }
}

median(spp_slp$slope, na.rm = T) #3.276331
mean(spp_slp$slope, na.rm = T) #4.406009
# NOTA!: Estos valores no distinguen a especies que ocupan la misma celda, osea, no hay un valor unico de pendiente por pixel, sino varios. De manera, que la mediana de la pendiente se debe sacar del raster de riqueza

# Pendiente por especie facet ---------------------------------------------

#unir la base de datos que contiene los valores de pendiente por pixel para cada especie, con los nombres de ordenes y familia en el objeto spp.nom.cods
spp_slp.join<-spp_slp |> 
  as_tibble() |> 
  merge(y= spp.nom.cods,
        by= "spp") |> 
  mutate(name=str_replace(name,"_"," "))

#write.table(spp_slp.join,"./outputs/tablas/spp_slp.join.txt", sep="\t", dec=".", row.names=F)

spp_slp.join<-read.table("./outputs/tablas/spp_slp.join.txt", dec=".", sep="\t", header=T)
head(spp_slp.join)

spp_slp.join |>
  # sample_n(100) |> 
  ggplot(aes(x=fct_reorder(name, slope, .desc=F, .na_rm = TRUE),y=slope))+
  geom_boxplot() +
  geom_hline(yintercept= 2.318058,color="red",linetype = "dashed")+
  geom_hline(yintercept = c(0,2,10,20), color="black") +
  labs(y="Terrain slope")+
  scale_y_continuous(trans="pseudo_log", limits=c(0, 35)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank())

ggsave(filename = "./outputs/F4A.png",
       width = 15,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 400)

# Cuantos pixeles en C y en T tiene cada una de las especies
spp_slp.join |> 
#  sample_n(100) |> 
  group_by(name, cat) |> 
  summarise(n=n()) |> 
  pivot_wider(names_from = cat, values_from = n) |> 
  mutate(dif=C-T) |> 
  arrange(dif) |> 
  pivot_longer(!c(name, dif),names_to = "cat", values_to = "count") |> 
  ggplot(aes(x=name, y=count, fill=cat))+
    geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank())
  



# cuantas especies se distribuyen debajo de la mediana de pendiente?

#mediana de la pendiente de mexico
mx.slope |> 
  as.data.frame() |> 
  summarise(median=median(slope_mx_g), # median: 1.585008
            mean=mean(slope_mx_g),
            sd=sd(slope_mx_g)) 

# Mediana de la pendiente potencialmente ocupada
Z.slp.cats |> 
  summarise(median=median(slope, na.rm=T), # median: 2.318058
            mean=mean(slope, na.rm=T),
            sd=sd(slope, na.rm=T)) 

# cuál es el valor de pendiente de cada especie en cada categoría?
spp_slp.join |> 
 #sample_n(size=100) |> 
  group_by(name, cat) |> 
  #  filter(endemismo=="Q") |> 
  summarise(median=median(slope, na.rm=T),
            sd=sd(slope, na.rm=T)) |> 
  pivot_wider(names_from = cat,
              values_from = c(median, sd)) |> 
  mutate(dif=median_C-median_T) |> 
  filter(dif<0) |> 
  arrange(dif) |> 
  print(n=60)

# cuál es la mediana de la pendiente en general para cada especie?
spp_slp.join |> 
  summarise(median=median(slope, na.rm=T), .by=name) |> 
  filter(median>2 & median<10) |> 
  arrange(median)

# 16 especies se distribuyeron en zonas planas (>2)
  16*100/56
# 40 especies se distribuyeron en zonas de ligera pendiente (2-10)
  40*100/56

# Cuáles especies sestán por encima del valor promedio de la pendiente?
  
  spp_slp.join |> 
    #  sample_n(size=100) |> 
    group_by(name) |> 
#    filter(endemismo=="Q") |> 
    summarise(median=median(slope, na.rm=T),
              sd=sd(slope, na.rm=T)) |> 
    filter(median<2.32) |> 
    arrange(median) |> 
    # filter(name=="Catharus olivascens" | 
    #        name=="Cardellina melanauris"|
    #       name=="Baeolophus wollweberi" |
    #       name=="Aphelocoma woodhouseii" |
    #       name=="Sitta pygmaea" |
    #       name=="Glaucidium hoskinsii") |> 
    print(n=60)
  
c("Catharus olivascens", "Cardellina melanauris", "Baeolophus wollweberi", "Aphelocoma woodhouseii", "Sitta pygmaea", "Glaucidium hoskinsii")
  # Total de especies
# Pendiente de todo el país
42*100/56 #42 por encima 
14*100/56 #14 por debajo 

# Pendiente del área ocupada
38*100/56 #38 por encima 
18*100/56 #18 por debajo 

  # Endemicas y cuasiendémicas
# Pendiente de todo el país
1*100/18 #1 por debajo
17*100/18 #17 por encima

# Pendiente de todo el país
2*100/18 #2 por debajo
16*100/18 #16 por encima


# Riqueza facet ordenado por slope ----------------------------------------
spp_slp.join |> 
 # sample_n(size=1000) |> 
  #  group_by(nom, slope) |> 
  #   summarise(count=median(slope)) |> 
  ggplot(aes(x=fct_reorder(str_replace(name,"_", " "), slope, .desc=F, .na_rm=T), y=slope, fill=cat)) +
  geom_boxplot() +
  geom_hline(aes(yintercept= mx.slope |>
                   as.data.frame() |>
                   as_tibble() |>
                   pull() |>
                   median()),color="red") +
  labs(y="Terrain slope", fill="Categories") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank()) +
  scale_fill_manual(values=c("#248f5d","#e5b636"), labels = c("Conserved","Transformed"))

# ggsave(filename = "./outputs/slp_spp_bxplt_facet.png",
#        width = 15,
#        height = 10, #alto
#        scale=2,
#        units ="cm",
#        dpi = 200)


# Area en cada valor de pendiente -----------------------------------------
mx.slope
area<-rast("./WorldClim_30s/wc2.1_30s_elev/Area_slp_cor_AG/AC_buf.tif")
INEGI_VII<-rast("./INEGI/Cambios/INEGI_VII_TC.tif")


#1. reproyectar capa de area
areawgs<-project(area, mx.slope, method="bilinear")
INEGI_VII<-resample(INEGI_VII, mx.slope, method="near")

#2. extraer valores de pendiente y area
slp.area<-terra::extract(c(areawgs,INEGI_VII), as.points(mx.slope), bind = T) |> 
  as.data.frame()

slp.area |> 
  na.omit() |> 
#  sample_n(size=100) |> 
  mutate(AC_buf=AC_buf/1000000, #a km cuadrado
         slope_clas=cut(slope_mx_g , breaks = 0:31, labels = 0:30)) |> 
  group_by(slope_clas, INEGI_VII_TC) |> 
  summarise(sum.area=sum(AC_buf)) |>
  na.omit() |> 
  #ggplot(aes(x=slope_clas, y=sum.area))+
   ggplot(aes(x=slope_clas, y=sum.area, fill=as.character(INEGI_VII_TC)))+
  geom_bar(stat='identity')+
  #labs(x="Terrain slope", y="Total area (km\u00b2)") +
   labs(x="Terrain slope", y="Total area (km\u00b2)", fill="Categories") +
   scale_fill_manual(values=c("#248f5d","#e5b636"), labels = c("Conserved","Transformed"))+
  theme_classic()

# ggsave(filename = "./outputs/slp_area2.png",
#        width = 15,
#        height = 10, #alto
#        scale=2,
#        units ="cm",
#        dpi = 200)


# Mapas bivariados --------------------------------------------------------

# Cargo la tabla para optimizar tiempo
#Para optimizar el tiempo de renderización del mapa, escalaré todo a 0.1 grados
#spp_PAM contiene las coordenadas de los centroides de los pixeles de distribución de cada especie
spp_PAM<- read.table("./outputs/tablas/spp_PAM_C.txt",sep="\t", dec=".", header=T) #OJO esto es con la riqueza en C
head(spp_PAM)

#Para poder hacer los modelos de regresión, es necesario crear un mapa de riqueza y con letsR puedo crear un mapa de riqueza apartir de una tabla de las coordanas de las especies (spd.list.join) 

#1. Hay que crear la Matriz de presencia-ausencia PAM a partir de los puntos de presencia a una resolucion de 10km
library(letsR)
PAM<-lets.presab.points(xy= as.matrix(spp_PAM[,1:2]),
                        species= spp_PAM$spp,
                        xmn = -120,
                        xmx = -85,
                        ymn = 14,
                        ymx = 33, #10km 0.08333 / 30km 0.24999 /50km 0.41665
                        resol = 0.1) 

x11()
plot(PAM)
summary(PAM)
PAM.raster<-PAM$Richness_Raster

#Hay un problema con letsR, me corta la 


# 2. Cargar y remuestrear el rasters a la resolución de la capa de riqueza ----------

# Remuestrear pendiente ---------------------------------------------------

mx.slope
slope.3<-resample(mx.slope, PAM.raster, method="mode")

#pendiente vs riqueza
slp.x.Z<-slope.3 |> 
  as.polygons(aggregate=F) |> 
  terra::extract(x=PAM.raster, bind=T) |> 
  st_as_sf()


slp.x.Z<-slp.x.Z[slp.x.Z$lyr.1!=0,]

# #intervalos de pendiente
# x=classInt::classIntervals(slp.x.Z$slope_mx_g,n=6,style = "fixed",
#                            fixedBreaks=c(0, 7, 13.9, 20.9,27.8,34.7))
# #intervalos de riqueza
# y=classInt::classIntervals(slp.x.Z$lyr.1,6,style = "fixed",
#                            fixedBreaks=c(0, 5, 10, 15, 20,30,40))
#intervalos de pendiente
x=classInt::classIntervals(slp.x.Z$slope_mx_g,n=4,style = "fixed",
                           fixedBreaks=c(0, 2, 10, 20, 34.7))
x
#intervalos de riqueza 
y=classInt::classIntervals(slp.x.Z$lyr.1,3,style = "fixed", intervalClosure = "right",
                           fixedBreaks=c(1, 3, 5, 10, 40))

y

slp.x.Z$x = classInt::findCols(x)
slp.x.Z$y = classInt::findCols(y)
slp.x.Z$alpha = as.character(slp.x.Z$x + slp.x.Z$y)
slp.x.Z$color = as.character(atan(slp.x.Z$y/slp.x.Z$x))

map<-
  ggplot()+
  geom_sf(data = slp.x.Z,aes(fill=color, alpha=alpha),shape=15, size=11,show.legend = F, lwd=0)+
  scale_fill_viridis_d(option="inferno")+
  theme_void()
map

  
#d=expand.grid(x=1:3,y=1:3)

leg<- 
  ggplot(expand.grid(x=1:4,y=1:4), aes(x,y))+
  geom_tile(aes(alpha=x+y,fill=atan(y/x)))+
  scale_fill_viridis_c(option="inferno")+
  labs(x="Terrain slope", y="Richness")+
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())+
  coord_equal()
leg

# library(patchwork)
map + inset_element(leg, 0,0,0.4,0.4)

ggsave(filename = "./outputs/vibariate.png",
       width = 15,
       height = 10, #alto
       scale=1,
       units ="cm",
       dpi = 200)

# Distribución de frecuencia de las combinaciones de riqueza y pendiente

slp.x.Z |> 
  as.data.frame() |> 
  rename(slope=slope_mx_g, z=lyr.1) |>
  mutate(xy=paste0(x,y)) |> 
  select(!geometry) |> 
  filter(z!=0) |> 
  summarise(freq=n(), .by=c(x,y,xy, color, alpha)) |>
  arrange(y) |> 
  mutate(color=seq(1:15)) |> 
  ggplot(aes(x=fct_reorder(xy, freq, .desc=T), y=freq, fill=as.character(xy)))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c("#fae9eb", "#ffe6bf", "#ffefc5", "#f9fcca","#d4c1d8", "#f1b0b9", "#ffb38a", "#ffbd4e","#a79da8", "#c080ac", "#e57686", "#ff7852","#737375", "#ad407f", "#d43852", "#875396")) +
  labs(y="Frecuencia de pixeles")+
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave(filename = "./outputs/vibariate_freq.svg",
       width = 12,
       height = 15, #alto
       scale=1,
       units ="cm",
       dpi = 200)

# lo mismo, pero en un grafico treemap
#install.packages("treemapify")
library(treemapify)

slp.x.Z |> 
  as.data.frame() |> 
  rename(slope=slope_mx_g, z=lyr.1) |>
  mutate(xy=paste0(x,y)) |> 
  select(!geometry) |> 
  filter(z!=0) |> 
  summarise(freq=n(), .by=c(x,y,xy, color, alpha)) |>
  mutate(fr.prc=freq*100/nrow(slp.x.Z),
         xy2=recode(xy, '11'='Po-Pl','12'='B-Pl','13'='I-Pl','14'='R-Pl','21'='Po-L','22'='B-L','23'='I-L','24'='R-L','31'='Po-M','32'='B-M','33'='I-M','34'='R-M','41'='Po-A','43'='I-A','44'='R-A')) |> 
  arrange(y) |> 
  ggplot(aes(area = freq, fill = xy, label = paste0(xy2,'\n',round(fr.prc,1),'%'))) +
  geom_treemap() +
  geom_treemap_text(place = "centre",
                    size = 12)+
  scale_fill_manual(values=c("#fae9eb", "#ffe6bf", "#ffefc5", "#f9fcca","#d4c1d8", "#f1b0b9", "#ffb38a", "#ffbd4e","#a79da8", "#c080ac", "#e57686", "#ff7852","#737375", "#ad407f", "#d43852", "#875396"))

ggsave(filename = "./outputs/vibariate_freq2.svg",
       width = 12,
       height = 15, #alto
       scale=1,
       units ="cm",
       dpi = 200)

#Fin del script
