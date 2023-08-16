#En este script voy a hacer todo lo que lleva el postprocesameinto de los modelos

# Binarización de los modelos ---------------------------------------------
Mods.dir<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = T) #modelos direcciones 
spp.mod.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F)

tabla.var.bin<-data.frame(spp=list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/M_spec.vbles_DD/", full.names = F), `10ptp.median`=NA, `10ptp.mean`=NA, `10ptp.max`=NA, `10ptp.rep.max`=NA, `10ptp.min`=NA, `10ptp.rep.min`
                          =NA, Mean_AUC_ratio=NA, pval_pROC=NA, Omiss.rat=NA, AICc=NA, delta_AICc=NA, W_AICc=NA, num_par=NA)


#max.res[,39] #X10.percentile.training.presence.Logistic.threshold 

#Primera parte: Extracción de 10 percentil training presence por especie y algunas métricas
#comparador
tab.mods<-read.table("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.lista.txt", sep="\t", dec=".", header=T)

#library(stringr)
for (x in 1:length(tab.mods[,1])){
  #1. arbir carpeta con modelos finales
  carp.fm<- list.files(paste0(Mods.dir[x], "/Final_Models"), full.names=F)
  #1.1 Aplicar if else porque el órden de las carpetas es diferente al nombre de las variables
  if (tab.mods[x,1]==substr(spp.mod.list,1,nchar(spp.mod.list)-2)[x]) { 
  #2. Calculo de la mediana
  #2.1 abrir maxent results
  max.res<-read.csv(paste0(paste0(Mods.dir[x], "/Final_Models/", tab.mods[x,2], "/maxentResults.csv")))
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
  cal.res<-read.csv(paste0(paste0(Mods.dir[x], "/Calibration_results/calibration_results.csv")))
  tabla.var.bin[x,8]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,2]
  tabla.var.bin[x,9]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,3]
  tabla.var.bin[x,10]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,4]
  tabla.var.bin[x,11]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,5]
  tabla.var.bin[x,12]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,6]
  tabla.var.bin[x,13]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,7]
  tabla.var.bin[x,14]<-cal.res[cal.res[,1] %in% tab.mods[x,2],][,8]
  }
}

#write.csv2(tabla.var.bin, "D:/Maestria_DD/spp.records_DD/Info.FM.csv", row.names = F, dec=".")
tabla.var.bin<-read.csv2("D:/Maestria_DD/spp.records_DD/Info.FM.csv")

#Segunda parte: Binarizar los modelos usando el umbral de cada especie
library(terra)

for (x in 1:8) { #length(tab.mods[,1])
#1. Abrir el raster por especie
ras.median<-rast(paste0(Mods.dir[x], "/Final_Models/", tab.mods[x,2], "/",
  paste(substr(spp.mod.list, 1,1), #obtengo la primera letra
        substr(word(gsub("_", " ", spp.mod.list),start=2,end=2,sep=fixed(" ")),1,4), sep="_")[x], 
  "_median.asc"))
ras.median2<-ras.median
ras.median2[ras.median2>tabla.var.bin[x,2]]<-1
ras.median2[ras.median2<=tabla.var.bin[x,2]]<-0
writeRaster(ras.median2, paste0("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", 
                                spp.mod.list[x] |> #nombre del .tif de cada especie
                                  word(1, sep = fixed("_")) |> #separo la primera palabra (genero)
                                  str_sub(1,2) |> #solo extraigo la primera letra
                                  paste( #creo el paste para unir genero y especie
                                    spp.mod.list[x] |> #nombre del .tif de cada especie
                                      word(2, sep = fixed("_")) |> #extraigo la segunda palabra
                                      substr(1,4), sep="_")#extraigo los 4 primeros caracteres
                                , "_bm.tif"), overwrite=TRUE) #binarian model
print(spp.mod.list[x])}

# -------------------------------------------------Mapas de riqueza --------------------------------------------------------

# Raster de riqueza total -------------------------------------------------

spp.bin.dir<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", pattern=".tif$", full.names = T) #ubicación de la carpeta con los Modelos_binarios
spp.bin.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", pattern=".tif$", full.names = F)
ext.ref<-rast("D:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/wc2.1_30s_bio_1_Ame.tif") #raster con la extension que quiero

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
  writeRaster("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_total.tif", 
              overwrite=T)

# Riqueza en México -------------------------------------------------------

#1. Recortar modelos a méxico
ext.ref.mx<-vect("D:/Maestria_DD/Shapes_MSc_DD/Mexico_Estados/Mexico_continent.shp")
spp.bin.dir #direccion de los Modelos_binarios
spp.bin.list #lista de nombres de modelos binarios
trans.vec<-vect("D:/Maestria_DD/Shapes_MSc_DD/INEGI/INEGI_II_702825007021_s/conjunto_de_datos/INEGI_II_ZU-DV.shp") #vector que contiene las zonas transformadas

for (x in 1:length(spp.bin.dir)) {
  spp.bin.dir[x] |>
    rast() |>
    mask(ext.ref, touches=F) |>
    crop(ext.ref) |> 
    mask(trans.vec, inverse=T, touches=F) |>  #corrección por las zonas transformadas
    writeRaster(
      paste0("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/",
        substr(spp.bin.list[x], 1, nchar(spp.bin.list[x]) - 4),
        "_MX.tif"),overwrite=T)
print(spp.bin.list[x])
}

#2. Calculo de riqueza
bin.list.mx<-list()
spp.bin.dir.mx<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = T) #ubicación de la carpeta con los Modelos_binarios
spp.bin.list.mx<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", pattern=".tif$", full.names = F)

# abro todos los rasters en una lista
for (i in 1:length(spp.bin.dir.mx)) {
  bin.list.mx[[i]]<-rast(spp.bin.dir.mx[i])
  names(bin.list.mx)[i]<-substr(spp.bin.list.mx[i], 1,nchar(spp.bin.list.mx[i])-4)
}

#calculo de riqueza
bin.list.mx |> 
  lapply(FUN=extend, y=ext(ext.ref.mx)) |> 
  rast() |> 
  sum(na.rm=T) |> 
  writeRaster("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif") |> 
  plot()

# ------------------------------------------------- TMC --------------------------------------------------------

# Corte a categorias ------------------------------------------------------

shp.cambios.dir<- list.files("D:/Maestria_DD/Shapes_MSc_DD/INEGI/Cambios/", pattern = ".shp$", full.names = T)
shp.cambios.dir<-shp.cambios.dir[-which(shp.cambios.dir=="D:/Maestria_DD/Shapes_MSc_DD/INEGI/Cambios/S2S7_cambio_shp.shp")]

cats2<-list.files("D:/Maestria_DD/Shapes_MSc_DD/INEGI/Cambios/", pattern = ".shp$", full.names = F)
cats<-substr(cats2, 1,nchar(cats2)-4)[-7]

library(stringr)
for (x in 1:length(shp.cambios.dir)) {# contador de las categorias
  cat.shp<-vect(shp.cambios.dir[x])
  cat.mask<-mask(ras.bin.inmask, cat.shp, touches=F) #touches=F solo los centroides que caigan dentro del poligono
  writeRaster(cat.mask, paste0("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TMC/Z_",cats[x],"cor.tif"), overwrite=T)
  }
#calcular numero de pixeles por categoría
cat.dir<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TMC/", pattern=".tif$", full.names = T)
cat.list<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TMC/", pattern=".tif$", full.names = F)

list.cat<-list()
for (i in 1:length(cat.dir)) {
  list.cat[[i]]<-rast(cat.dir[i])
  names(list.cat)[i]<-substr(cat.list[i], 1,nchar(cat.list[i])-4)
}

bl.stack.cat<-rast(lapply(list.cat, FUN=extend, y=ext(ext.ref)))
count<-freq(bl.stack.cat)

count2<-data.frame(tapply(count$count, count$layer, sum))
row.names(count2)<-names(bl.stack.cat)
names(count2)<-"count pix"


# Relación con la pendiente -----------------------------------------------
mx.slope<- rast("D:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_slope_MX.tif")

#1. convertir raster de categorías a puntos

cat.Z.list<-list()
for(x in 1:length(cat.dir)){
  #1. convertir raster de riqueza por categorías en puntos
  cat.point<-as.points(bl.stack.cat[[x]])
  #2. Extraer los valores de pendiente por valor de riqueza
  tabla<-as.data.frame(extract(mx.slope, cat.point, bind=T)) #contiene valores de cero
  tabla$cat<-substr(cat.list[x], 3,nchar(cat.list[x])-7)
  names(tabla)<-c("riqueza", "slope", "cat")
  #2.1. eliminar valores de cero
  tabla.subs<-subset(tabla, riqueza!=0)
  #2.2 guardarlo en una lista
  cat.Z.list[[x]]<-tabla.subs
  names(cat.Z.list)[x]<-paste0(substr(cat.list[x], 9,nchar(cat.list[x])-4),"_Z")
  #2.3 guardar el txt
  write.table(tabla.subs, paste0("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TMC/", levels(factor(tabla$cat)), "_Z_slp.txt"), sep="\t", dec=".", row.names = F)
}

# Plots de densidad -------------------------------------------------------
dir.z<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/TMC/", pattern=".txt$", full.names = T)
CC_11_Z<-read.table(dir.z[1], sep="\t", dec=".", header=T)
CM_12_Z<-read.table(dir.z[2], sep="\t", dec=".", header=T)
CT_13_Z<-read.table(dir.z[3], sep="\t", dec=".", header=T)
MC_21_Z<-read.table(dir.z[4], sep="\t", dec=".", header=T)
MM_22_Z<-read.table(dir.z[5], sep="\t", dec=".", header=T)
MT_23_Z<-read.table(dir.z[6], sep="\t", dec=".", header=T)
TC_31_Z<-read.table(dir.z[7], sep="\t", dec=".", header=T)
TM_32_Z<-read.table(dir.z[8], sep="\t", dec=".", header=T)
TT_33_Z<-read.table(dir.z[9], sep="\t", dec=".", header=T)

#unirlos todos en un solo DFZ
lapply(cat.Z.list, summary)
Z.table.cats<-as.data.frame(do.call(rbind, cat.Z.list))
Z.table.cats$count<-1
Z.table.cats<-Z.table.cats[-which(is.na(Z.table.cats[,2])),]
summary(Z.table.cats)

##Aqui voy a graficar el numero de pixeles por valor de pendiente
CC_11_Z$count<-1
summary(Z.table.cats)

tapply(Z.table.cats$slope, Z.table.cats$cat, sd)
row.names(count2)<-names(bl.stack.cat)
names(count2)<-"count pix"

#un grafico de violin
ggplot(data=Z.table.cats, aes(y=slope, x=cat, fill=cat))+
  geom_violin() +
  geom_boxplot(width=0.15) + theme_minimal()
  
#densidades de kernell
CC_11_Z.sna<-CC_11_Z[-which(is.na(CC_11_Z$slope)),]
CM_12_Z.sna<-CM_12_Z[-which(is.na(CM_12_Z$slope)),]
CT_13_Z.sna<-CT_13_Z[-which(is.na(CT_13_Z$slope)),]
MC_21_Z.sna<-MC_21_Z[-which(is.na(MC_21_Z$slope)),]
MM_22_Z.sna<-MM_22_Z[-which(is.na(MM_22_Z$slope)),]
MT_23_Z.sna<-MT_23_Z[-which(is.na(MT_23_Z$slope)),]
TC_31_Z.sna<-TC_31_Z[-which(is.na(TC_31_Z$slope)),]
TM_32_Z.sna<-TM_32_Z[-which(is.na(TM_32_Z$slope)),]
TT_33_Z.sna<-TT_33_Z[-which(is.na(TT_33_Z$slope)),]

plot(density(CC_11_Z.sna$slope), ylim=c(0,0.7), xlim=c(0,70))
lines(density(CM_12_Z.sna$slope), col="red")
lines(density(CT_13_Z.sna$slope), col="green")
lines(density(MC_21_Z.sna$slope), col="blue")
lines(density(MM_22_Z.sna$slope), col="yellow")
lines(density(MT_23_Z.sna$slope), col="gray")
lines(density(TC_31_Z.sna$slope), col="skyblue")
lines(density(TM_32_Z.sna$slope), col="slateblue")
lines(density(TT_33_Z.sna$slope), col="sienna")

#plot total
ggplot(Z.table.cats, aes(x=slope, fill=cat, color=cat)) + 
  geom_density(alpha=0.2, linewidth=1) +
  scale_fill_manual(values=c("black","red","green","blue","yellow","gray","skyblue","aquamarine4","sienna")) + 
  scale_color_manual(values=c("black","red","green","blue","yellow","gray","skyblue","aquamarine4","sienna"))

#Plot por cada uno
CC_11.denp<-ggplot(CC_11_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="black") + 
  scale_color_manual(values="black") +
  labs(title="Conservado - Conservado") +
  theme(legend.position="none") +
  ylim(0,0.3)

CM_12.denp<-ggplot(CM_12_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="red") + 
  scale_color_manual(values="red")+
  labs(title="Conservado - Medio") +
  theme(legend.position="none")+
  ylim(0,0.3)

CT_13.denp<-ggplot(CT_13_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="green") + 
  scale_color_manual(values="green")+
  labs(title="Conservado - Transformado") +
  theme(legend.position="none")+
  ylim(0,0.3)

MC_21.denp<-ggplot(MC_21_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="blue") + 
  scale_color_manual(values="blue")+
  labs(title="Medio - Conservado") +
  theme(legend.position="none")+
  ylim(0,0.3)

MM_22.denp<-ggplot(MM_22_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="yellow") + 
  scale_color_manual(values="yellow")+
  labs(title="Medio - Medio") +
  theme(legend.position="none")+
  ylim(0,0.3)

MT_23.denp<-ggplot(MT_23_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="gray") + 
  scale_color_manual(values="gray")+
  labs(title="Medio - Transformado") +
  theme(legend.position="none")+
  ylim(0,0.3)

TC_31.denp<-ggplot(TC_31_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="skyblue") + 
  scale_color_manual(values="skyblue") +
  labs(title="Transformado - Conservado") +
  theme(legend.position="none")+
  ylim(0,0.3)

TM_32.denp<-ggplot(TM_32_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="aquamarine4") + 
  scale_color_manual(values="aquamarine4") +
  labs(title="Transformado - Medio") +
  theme(legend.position="none")+
  ylim(0,0.3)

TT_33.denp<-ggplot(TT_33_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="sienna") + 
  scale_color_manual(values="sienna") +
  labs(title="Transformado - Transformado") +
  theme(legend.position="none")+
  ylim(0,0.3)

#library(gridExtra)
grid.arrange(CC_11.denp,
             CM_12.denp,
             CT_13.denp,
             MC_21.denp,
             MM_22.denp,
             MT_23.denp,
             TC_31.denp,
             TM_32.denp,
             TT_33.denp)


# Plot de riqueza por categoria -------------------------------------------
CC_11.denp<-ggplot(CC_11_Z.sna, aes(x=riqueza, color=cat, fill=cat)) + 
  geom_density(alpha=0.5, linewidth=0.9) +
  scale_fill_manual(values="black") + 
  scale_color_manual(values="black") +
  labs(title="Conservado - Conservado") +
  theme(legend.position="none") +
  ylim(0,0.3)

# Plots de Riqueza vs pendiente--------------------------------------------------------
lapply(list(CC_11_Z.sna,
CM_12_Z.sna,
CT_13_Z.sna,
MC_21_Z.sna,
MM_22_Z.sna,
MT_23_Z.sna,
TC_31_Z.sna,
TM_32_Z.sna,
TT_33_Z.sna), FUN= head)

#plot x ~ y
plot(CC_11_Z.sna$slope, CC_11_Z.sna$riqueza)
a<-lm(CC_11_Z.sna$riqueza ~ CC_11_Z.sna$slope)
abline(a, color="red")

CC_11_Z.grap<- ggplot(data=CC_11_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(CC_11_Z.sna$slope), color="red") +
  labs(title="Conservado - Conservado") + 
  ylim(0,69.04134)

CM_12_Z.grap<- ggplot(data=CM_12_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(CM_12_Z.sna$slope), color="red") +
  labs(title="Conservado - Medio") + 
  ylim(0,69.04134)

CT_13_Z.grap<- ggplot(data=CT_13_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(CT_13_Z.sna$slope), color="red") +
  labs(title="Conservado - Transformado") + 
  ylim(0,69.04134)

MC_21_Z.grap<- ggplot(data=MC_21_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(MC_21_Z.sna$slope), color="red") +
  labs(title="Medio - Conservado") + 
  ylim(0,69.04134)

MM_22_Z.grap<- ggplot(data=MM_22_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(MM_22_Z.sna$slope), color="red") +
  labs(title="Medio - Medio") + 
  ylim(0,69.04134)

MT_23_Z.grap<- ggplot(data=MT_23_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(MT_23_Z.sna$slope), color="red") +
  labs(title="Medio - Transformado") + 
  ylim(0,69.04134)

TC_31_Z.grap<- ggplot(data=TC_31_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(TC_31_Z.sna$slope), color="red") +
  labs(title="Transformado - Conservado") + 
  ylim(0,69.04134)

TM_32_Z.grap<- ggplot(data=TM_32_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(TM_32_Z.sna$slope), color="red") +
  labs(title="Transformado - Medio") + 
  ylim(0,69.04134)

TT_33_Z.grap<- ggplot(data=TT_33_Z.sna, aes(x=factor(riqueza), y=slope))+
  geom_boxplot() +
  geom_hline(yintercept = median(TT_33_Z.sna$slope), color="red") +
  labs(title="Transformado - Transformado") + 
  ylim(0,69.04134)

grid.arrange(CC_11_Z.grap,
             CM_12_Z.grap,
             CT_13_Z.grap,
             MC_21_Z.grap,
             MM_22_Z.grap,
             MT_23_Z.grap,
             TC_31_Z.grap,
             TM_32_Z.grap,
             TT_33_Z.grap)

geom_text(aes(x=slope, fill=cat,label = cat))
