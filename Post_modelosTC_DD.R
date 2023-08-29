#En este script voy a hacer todo lo que lleva el postprocesamiento de los modelos   

#-------------------------------------------------------------------------------------------------------------- # 

# ------------------------------------- Binarización de los modelos ---------------------------------------------

#-------------------------------------------------------------------------------------------------------------- # 

mods.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = T) #modelos direcciones 
spp.mod.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F)

tabla.var.bin<-data.frame(spp=list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F), `10ptp.median`=NA, `10ptp.mean`=NA, `10ptp.max`=NA, `10ptp.rep.max`=NA, `10ptp.min`=NA, `10ptp.rep.min`
                          =NA, Mean_AUC_ratio=NA, pval_pROC=NA, Omiss.rat=NA, AICc=NA, delta_AICc=NA, W_AICc=NA, num_par=NA)

#max.res[,39] #X10.percentile.training.presence.Logistic.threshold 

#Primera parte: Extracción de 10 percentil training presence por especie y algunas métricas
#comparador
tab.mods<-read.table("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.lista.txt", sep="\t", dec=".", header=T)

#library(stringr)
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
library(terra)

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

#-------------------------------------------------------------------------------------------------------------- # 

# ----------------------------------------------- T y C S7 ----------------------------------------------------  

#-------------------------------------------------------------------------------------------------------------- #

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
cat.dir.TC<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", full.names = T)[-c(2,4)] # dirección de las carpetas de categorias
cat.list.TC<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/",full.names = F)[-c(2,4)] #lista de categorias
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

sum(tab.distxcats[1,12:13]) #100% ok

tab.distxcats[which(tab.distxcats$endemismo=="E"),3]<-1
tab.distxcats[which(tab.distxcats$endemismo=="Q"),3]<-2
tab.distxcats[which(tab.distxcats$endemismo=="NE"),3]<-3

color<- tab.distxcats$value
  
color[which(color==1)]<-"lightblue"
color[which(color==2)]<-"#A01127"
color[which(color==3)]<-"#4D4D4D"


tab.distxcats |> 
  pivot_longer(cols=c(pc.dist.C, pc.dist.T),
               names_to= "pc.cats",
               values_to = "pc") |> 
   # group_by(nom, pc.cats) |> 
   # summarise(count=sum(pc)) |> 
  ggplot(aes(x=fct_reorder(str_replace(name,"_", " "), value,.desc=F), y=pc, fill=pc.cats))+
  geom_bar(stat="identity") +
  labs(y="Proporcion", fill="Categorias") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 13),
        axis.text.y = element_text(size = 13),
        panel.background = element_blank()) +
  geom_hline(yintercept = 50, col="black", linetype="dashed")+
  scale_fill_manual(values=c("#248f5d","#f56038"), labels = c("C","T"))

ggsave(file = "./outputs/catsxspp.svg",
       width = 1706,
       height = 1271,
       scale=3,
       units ="px")

#---------------------------------------- Relación con la pendiente -----------------------------------------------
mx.slope<- rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
hfp<-rast("F:/Maestria_DD/Shapes_MSc_DD/HFP/Dryad_2020/Human_footprint_maps/hfp2013_merisINT_MX.tif")
ecoreg<- rast("F:/Maestria_DD/Shapes_MSc_DD/Provincias biogeograficas/Ecoregiones/ecoreg.tif")

#El extent no es el mismo!!!!
# resample(mx.slope, mx.slope2, method="bilinear")|> 
#   writeRaster("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif", overwrite=T)

TC.Z.dir<- list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", pattern=".tif$", full.names = T) # Direcciones de la riqueza por categoría
TC.Z.list<- list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", pattern=".tif$", full.names = F) # Lista de los nombres de los rasters de riqueza por categoría
cat.dir.TC

bl.stack.cat.TC<- rast(TC.Z.dir)
#1. convertir raster de categorías a puntos

cat.Z.list.TC<-list()
for(x in 1:length(cat.dir.TC)){
  #1. convertir raster de riqueza por categorías en puntos
  cat.point<-as.points(bl.stack.cat.TC[[x]]) #contiene los rasters de riqueza por categoria
  #2. Extraer los valores de pendiente por valor de riqueza
  c.p<-terra::extract(mx.slope, cat.point, bind=T, xy=T)  
  c.p.hfp<-terra::extract(x=hfp, y=c.p, bind=T) #tambien el hfp 
  #2.1 Convertir tabla en un DF
  tabla<-as.data.frame(terra::extract(x=ecoreg, y=c.p.hfp, bind=T)) #contiene valores de cero
  tabla$cat<-substr(cat.list.TC[x], 3,4)
  names(tabla)<-c("z", "slope", "x","y", "hfp", "ecoreg", "cat")
  tabla<-tabla[c(3,4,1,2,5,6,7)]
  tabla.subs<-subset(tabla, z!=0)
  #2.2 guardarlo en una lista
  cat.Z.list.TC[[x]]<-tabla.subs
  names(cat.Z.list.TC)[x]<-paste0(substr(cat.list.TC[x], 3,4),"_Z")
  #2.3 guardar el txt
  write.table(tabla.subs, paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TC_s7/", levels(factor(tabla$cat)), "_Z_slp.txt"), sep="\t", dec=".", row.names = F)
}

#unirlos todos en un solo DF
Z.slp.cats<-as.data.frame(do.call(rbind, cat.Z.list.TC)) #contiene todos los valores de pendiente por valor de riqueza en cada categoria
Z.slp.cats$count<-1 #agrego una columan de contador
Z.slp.cats<-Z.slp.cats |> 
  #elimino NA
  na.omit() # |> write.table("./outputs/Z.slp.cats.txt", sep = "\t", dec=".", row.names = F)
summary(Z.slp.cats) #verifico

Z.slp.cats<- read.table("./outputs/Z.slp.cats.txt", sep = "\t", dec=".", header=T)

# densidades de kernel ----------------------------------------------------

# Slope -------------------------------------------------------------------
#plot total
ggplot(Z.slp.cats, aes(x=slope, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#f56038"),labels=c("Conserved", "Transformed")) + 
  scale_color_manual(values=c("#248f5d","#f56038"), labels=c("Conserved", "Transformed")) +
  labs(x="Terrain slope", y="Density", fill="Categories", color="Categories")

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
  scale_fill_manual(values=c("#248f5d","#f56038")) + 
  scale_color_manual(values=c("#248f5d","#f56038")) +
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


# Riqueza -----------------------------------------------------------------
#plot total
ggplot(Z.slp.cats, aes(x=z, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#f56038"),labels=c("Conserved", "Transformed")) + 
  scale_color_manual(values=c("#248f5d","#f56038"), labels=c("Conserved", "Transformed")) +
  labs(x="Richness", y="Density", fill="Categories", color="Categories")

# ggsave(filename = "./outputs/riq_dens_tot.png",
#          width = 15,
#          height = 10, #alto
#          scale=2,
#          units ="cm",
#          dpi = 200)

#por cada uno
f_labels
#f_Z<- #esto es para el plot doble
ggplot(Z.slp.cats, aes(x=z, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("#248f5d","#f56038")) + 
  scale_color_manual(values=c("#248f5d","#f56038")) +
  facet_wrap("cat", labeller = labeller(cat=c(C= "Conserved",
                                              T= "Transformed"))) +
  theme_classic()+
  labs(x="Richness", y="Density")+
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

library(janitor)
ggplot(Z.slp.cats, aes(x=factor(z), y=slope)) + 
  geom_boxplot() +
  # group_by(nom, pc.cats) |> 
  # summarise(count=sum(pc)) |> 
  geom_hline(data= Z.slp.cats |> #esta linea crea un objeto con las medianas para cada facet
               group_by(cat) |> 
               summarize(mediana_slp = median(slope)),
             aes(yintercept= mediana_slp),
             color="red") +
  facet_wrap("cat", labeller = labeller(cat=c(C= "Conserved",
                                              T="Transformed")), nrow=2, ncol=1) +
  labs(x="Richness", y="Terrain slope") +
  theme_classic()+
  theme(legend.position="none",
        strip.text =  element_text(size=10, face="bold"))

ggsave(filename = "./outputs/slp_bxplt_facet.png",
       width = 10,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)


# barplot -----------------------------------------------------------------
Z.slp.cats |> 
  group_by(cat, z) |> 
  summarise(median= median(slope)) |> 
  ggplot(aes(x=z, y=median, color=cat))+
  geom_line() +
  facet_wrap("cat", labeller = labeller(cat=c(C= "Conserved",
                                              T= "Transformed")))+
  scale_color_manual(values=c("#248f5d","#f56038")) +
  labs(x="Richness", y="Terrain slope (median)", fill="Categories", color="Categories") +
  theme(legend.position = "none",
        strip.text =  element_text(size=10, face="bold"),
        strip.background = element_rect(color = "black", fill="white"))

#Las especies se concentran en valores más altos de pendiente

# ------------------------------------------- Pendente por especie ----------------------------------------

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

as.data.frame(do.call(rbind, list.cat.slp)) |> 
  write.table("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/spp_slp.txt", dec=".", sep="\t", row.names = F)
rm(list.cat.slp)

#Barplot
spp_slp<-read.table("F:/Maestria_DD/spp.records_DD/specialists_DD/tablas.outputs/spp_slp.txt", dec=".", sep="\t", header=T)
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

# Pendiente por especie facet ---------------------------------------------

#unir la base de datos que contiene los valores de pendiente por pixel para cada especie, con la información asociada en el objeto spp.nom.cods
spp_slp.join<-spp_slp |> 
  as_tibble() |> 
  merge(y= spp.nom.cods,
        by= "spp") 

spp_slp.join |>
  group_by(cat, name) |> 
  mutate(name=str_replace(name,"_"," ")) |> 
  summarize(median=median(slope, na.rm=T),
            sd=sd(slope, na.rm=T)) |> 
  ggplot(aes(x=reorder(name, median),y=median, fill=cat))+
  geom_bar(stat ="identity", position=position_dodge()) +
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank())+
  labs(y="Terrain slope (median)",fill="Categories")+
  scale_fill_manual(values=c("#248f5d","#f56038"), labels = c("Conserved","Transformed"))

  
# ggsave(filename = "./outputs/slp_spp_bxplt_facet2.png",
#        width = 15,
#        height = 10, #alto
#        scale=2,
#        units ="cm",
#        dpi = 200)

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
  scale_fill_manual(values=c("#248f5d","#f56038"), labels = c("Conserved","Transformed"))

ggsave(filename = "./outputs/slp_spp_bxplt_facet.png",
       width = 15,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)


# Mapas bivariados --------------------------------------------------------
#esto esta en veremos porque pienso incluir los graficos range-diversity de soberon (2021)
mx.slope
Riqueza_MX_HFP<-rast("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX_HFP.tif")
library(biscale)
library(sf)

#pendiente vs riqueza
slp.x.Z<- mx.slope |> 
  as.points() |> 
  extract(x=Riqueza_MX_HFP, bind=T) |> 
  st_as_sf()

slp.x.Z

head(slp.x.Z)

x=classInt::classIntervals(slp.x.Z$wc2.1_30s_slope,4,style = "quantile")
y=classInt::classIntervals(slp.x.Z$sum,4,style = "quantile")

slp.x.Z$x = classInt::findCols(x)
slp.x.Z$y = classInt::findCols(y)
slp.x.Z$alpha = as.character(slp.x.Z$x + slp.x.Z$y)
slp.x.Z$color = as.character(atan(slp.x.Z$y/slp.x.Z$x))

ggplot()+
  geom_sf(data = slp.x.Z,aes(fill=color,alpha=alpha),shape=15, size=11,show.legend = FALSE)+
  scale_fill_viridis_d()+
  theme_void()

leg<- ggplot(d, aes(x,y))+
  geom_tile(aes(alpha=x+y,fill=atan(y/x)))+
  scale_fill_viridis_c()+
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  coord_equal()

head(slp.x.Z)
x

d=expand.grid(x=1:3,y=1:3)
d

#Fin del script

