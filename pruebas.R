# Mapas bivariados --------------------------------------------------------

# Cargo la tabla para optimizar tiempo
#Para optimizar el tiempo de renderización del mapa, escalaré todo a 0.1 grados
#spp_PAM contiene las coordenadas de los centroides de los pixeles de distribución de cada especie
spp_PAM<- read.table("./outputs/tablas/spp_PAM_C.txt",sep="\t", dec=".", header=T)
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
mx.slope<- rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
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

### AQUI EMPIEZA LA DIFERENCIA

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

####### GRAFICOS DE MOSAICO

library(ggmosaic)
flights <- 
  fly  %>%
  filter(!is.na(do_you_recline), !is.na(rude_to_recline))

flights |> 
  select(rude_to_recline)

slp.x.Z |> 
  as.data.frame() |> 
  rename(slope=slope_mx_g, z=lyr.1) |>
  mutate(xy=paste0(x,y)) |> 
  select(!geometry) |> 
  filter(z!=0) |> 
  ggplot() +
  geom_mosaic(aes(x = product(x)), divider = "vspine")+ #x="slope", y="Riqueza"
  facet_grid(~y)+
  theme(aspect.ratio = 3,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

slp.x.Z |> 
  as.data.frame() |> 
  rename(slope=slope_mx_g, z=lyr.1) |>
  mutate(xy=paste0(x,y)) |> 
  select(!geometry) |> 
  filter(z!=0) |> 
  summarise(freq=n(), .by=c(x,y,xy, color, alpha)) |>
  mutate(fr.prc=freq*100/nrow(slp.x.Z)) |> 
  arrange(y) |> 
  mutate(color=seq(1:15)) |> 
  ggplot(aes(x=x,y=y)) +
  geom_point(aes(size=freq)) +
  scale_size_continuous(range = c(1, 20))


 # labs(x="Pendiente del terreno", y="Riqueza")

## donut chart
slp.x.Z |> 
  as.data.frame() |> 
  rename(slope=slope_mx_g, z=lyr.1) |>
  mutate(xy=paste0(x,y)) |> 
  select(!geometry) |> 
  filter(z!=0) |> 
  summarise(freq=n(), .by=c(x,y,xy, color, alpha)) |>
  mutate(fr.prc=freq*100/nrow(slp.x.Z),
         hsize=1) |> 
  arrange(y) |> 
  ggplot(aes(x = hsize, y = fr.prc, fill = xy)) +
  geom_col() +
  coord_polar(theta = "y")+
  xlim(c(0.2, 1 + 0.5)) +
  scale_fill_manual(values=c("#fae9eb", "#ffe6bf", "#ffefc5", "#f9fcca","#d4c1d8", "#f1b0b9", "#ffb38a", "#ffbd4e","#a79da8", "#c080ac", "#e57686", "#ff7852","#737375", "#ad407f", "#d43852", "#875396")) +
  geom_text(aes(label = fr.prc),
            position = position_stack(vjust = 0.5))+
  theme_classic()
  
install.packages("treemapify")
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

ggsave(filename = "./outputs/vibariate_freq.svg",
       width = 12,
       height = 15, #alto
       scale=1,
       units ="cm",
       dpi = 200)

# Zonificación de la pendiente en T y C ------------------------------------


library(terra)
library(tidyverse)
reg.bio<-vect("Provincias biogeograficas/Regiones_biog/Regiones_biog.shp")
biomas<-vect("Provincias biogeograficas/Ecoregions2017_dinerstein/Biomas_2017_Dinerstein.shp")
cat.rast<-rast("INEGI/Cambios/INEGI_VII_TC.tif")
z.mx<-rast("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif")
mx.slope<- rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
plot(cat.rast)

#esto puede ser para las regiones o para los biomas dependiendo de que quiera sacar, para no repetir el código que es el mismo
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

zon<-read.table("./outputs/tablas/zon.biom.txt", sep="\t", dec=".", header=T) #ojo debe cambiar el ojeto dependiendo de la zonificacion

zon |> 
  sample_n(size=200) |> 
  filter(z!=0) |> 
  ggplot(aes(x=cat, y=slope, fill=cat))+
  geom_boxplot()+
  scale_fill_manual(values=c("#248f5d","#e5b636"))+
  scale_y_continuous(trans="pseudo_log", limits=c(0, 35))+
  facet_wrap(~zonif)+
  geom_hline(yintercept= 2.318058,color="red",linetype = "dashed")+
  geom_hline(yintercept = c(0,2,10,20), color="black") +
  theme_classic()+
  labs(y="Pendiente del terreno")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())




#### ESTO ES DE OTRA VUELTA

#install.packages("virtualspecies")
library(virtualspecies)
#library(terra)

worldclim <- rast(list.files("F:/Maestria_DD/Proyectos/Nichos/Script/", pattern=".tif", full.names = T))

x1<-seq(1:100)
dnorm(1:300, mean = 250, sd=50) |> 
  plot()

betaFun(1:40, p1=17, p2 = 25, alpha = 0.6, gamma = 0.04) |> 
  plot(type="l") |> 
  abline(v=c(17,25), col="red")

my.parameters <- formatFunctions(wc2.1_10m_bio_1 = c(fun = 'betaFun', p1=17, p2 = 25, alpha = 0.9, gamma = 0.08),
                                 wc2.1_10m_bio_12 = c(fun = 'dnorm', mean = 4000, sd = 2000))

my.first.species <- generateSpFromFun(raster.stack = worldclim,
                                      parameters = my.parameters,
                                      plot = T)


my.first.species$suitab.raster |> 
  plot()

####

library(Rlab)

#dbern

# x values for the dbern() function
x <- c(0, 1, 0, 5, 7, 10)

# Using dbern() function to obtain the corresponding Bernoulli PDF
y <- dbern(x, prob = 0.5)

# Plotting dbern values
plot(x, y, type = "o")


# pbern -------------------------------------------------------------------

# x values for the
# pbern( ) function
x <- seq(0, 10, by = 1)

# using pbern( ) function
# to x to obtain corresponding
# Bernoulli  CDF
y <- pbern(x, prob = 0.5)  

# plot pbern values
plot(y, type = "o")    
