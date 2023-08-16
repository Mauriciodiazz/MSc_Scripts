###este sript contiene pruebas y analisis preliminares

### el compound topographic index es una función de la pendiente y de la contribución del flujo de agua por unidad de área, basicamente es una métrica de escorrentía
#Mi pregunta es: si yo hago una regresión o una  prueba estadistica de relación entre la pendiente y el cti ¿será significativa?

####especies virtuales????
library(virtualspecies)
library(raster)

vars<-stack(list.files("./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/", full.names=T)[1:2])
#supongamos que el rango de mi especie virtual, es el mimso rango de una especie real, pero, se va a distribuir de manera normal en todo su rango ambiental
library(terra)
a.wood<-vect("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Aphelocoma woodhouseii.shp")
a.ext<-extract(rast(vars[[1]]), a.wood)
head(a.ext)
max(a.ext$wc2.1_30s_bio_1)
hist(a.ext$wc2.1_30s_bio_1)
abline(v=mean(a.ext$wc2.1_30s_bio_1)) #mean=11.92878

x <- seq(-100, 100, by = .1)
y <- dnorm(x, mean = 11.92878, sd = 20)
plot(x,y)

shapiro.test(dnorm(x, mean = 12, sd = 20))

#definicion de parametros ambientales
my.parameters <- formatFunctions(wc2.1_30s_bio_1 = c(fun = 'dnorm', mean = 250, sd = 50),
                                 wc2.1_30s_bio_10 = c(fun = 'dnorm', mean = 4000, sd = 2000))


dnorm(x = 150, mean = 250, sd = 50)


a<-c(fun = 'dnorm', mean = 250, sd = 50)
a
spp.1


S2<-rast("D:/Maestria_DD/Shapes MSc_DD/INEGI/Cambios/INEGI_II_TC.tif")
S7<-rast("D:/Maestria_DD/Shapes MSc_DD/INEGI/Cambios/INEGI_VII_TC.tif")

s2s7<-sum(S2,S7)
plot(s2s7)
writeRaster(s2s7, "D:/Maestria_DD/Shapes MSc_DD/INEGI/Cambios/S2S7_CAMBIO_TC.tif",overwrite=TRUE)

library(skimr)
a<-read.table("Correlaciones/tabla_corr_es_pe.txt", sep="\t", dec=".", header=T)
skim(a)

library(tidyverse)
setwd("C:/Users/Mauricio Diaz/Documents/Maestria/")
files<-list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = T)
files

files.to<-
  files |> 
  str_replace_all(pattern=" ", replacement = "_")
files.to

file.rename(files, files.to)


a<-data.frame(nom.com=spp.mod.list, nom.cod= spp.mod.list |> #nombre del .tif de cada especie
                word(1, sep = fixed("_")) |> #separo la primera palabra (genero)
                str_sub(1, 2) |> #solo extraigo la primera letra
                paste(#creo el paste para unir genero y especie
                  spp.mod.list |> #nombre del .tif de cada especie
                    word(2, sep = fixed("_")) |> #extraigo la segunda palabra
                    substr(1, 4), sep = "_"))

list.files("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/", pattern = ".tif$") |> 
  data.frame()




a<-data.frame(nom=c("Aphelocoma_woodhouseii", "Atlapetes_pileatus", "Attila_pacificus", "Baeolophus_wollweberi", "Basileuterus_lachrymosus", "Cardellina_melanauris", "Cardellina_rubra", "Caryothraustes_poliogaster", "Catharus_dryas", "Catharus_mexicanus", "Catharus_occidentalis", "Catharus_olivascens", "Chlorestes_candida", "Claravis_pretiosa", "Coccothraustes_abeillei", "Crypturellus_occidentalis", "Dendrocolaptes_sheffleri", "Dysithamnus_mentalis", "Eucometis_penicillata", "Glaucidium_hoskinsii", "Habia_fuscicauda", "Lampornis_amethystinus", "Lanio_aurantius", "Lepidocolaptes_affinis", "Lepidocolaptes_leucogaster", "Leptotila_verreauxi", "Lipaugus_unirufus", "Lophornis_helenae", "Loxia_stricklandi", "Malacoptila_panamensis", "Megascops_seductus", "Micrastur_semitorquatus", "Microcerculus_philomela", "Microrhopias_quixensis", "Mionectes_oleagineus", "Myadestes_occidentalis", "Notharchus_hyperrhynchus", "Oreophasis_derbianus", "Pachyramphus_cinnamomeus", "Penelope_purpurascens", "Picoides_stricklandi", "Pipilo_nigrescens", "Poecile_sclateri", "Pseudastur_albicollis", "Sclerurus_mexicanus", "Sitta_pygmaea", "Strix_fulvescens", "Thamnistes_anabatinus", "Troglodytes_rufociliatus", "Trogon_ambiguus", "Trogon_mexicanus", "Tunchiornis_ochraceiceps", "Uropsila_pacifica", "Vireo_hypochryseus", "Vireo_paluster", "Zentrygon_carrikeri"),
              nom.cods= a$nom|> 
                word(1, sep= fixed("_")) |> 
                str_sub(1,2) |> 
                paste(a$nom |> 
                        word(2, sep=fixed("_")) |> 
                        substr(1,4), sep="_"),
           end=c("NE", "E", "E", "NE", "NE", "E", "E", "NE", "NE", "NE", "E", "E", "NE", "NE", "Q", "E", "E", "NE", "NE", "E", "NE", "E", "NE", "NE", "E", "NE", "NE", "NE", "Q", "NE", "E", "NE", "NE", "NE", "NE", "NE", "NE", "Q", "NE", "NE", "E", "E", "Q", "NE", "NE", "NE", "NE", "NE", "NE", "Q", "NE", "NE", "E", "E", "E", "E"))

a <- rast("D:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_total.tif")
summary(a)
b<- vect("D:/Maestria_DD/Shapes_MSc_DD/HFP/Dryad_2020/Human_footprint_maps/hfp2013_merisINT_MX_20.shp")
mex<- vect("D:/Maestria_DD/Shapes_MSc_DD/Mexico_Estados/Mexico_continent.shp")


a |> 
  mask(mask = mex, touches=F) |> 
  mask(mask = b, touches=F, inverse=T) |> 
  writeRaster("./outputs/rastercorr.tif")





# Mapa bivariado ----------------------------------------------------------
#ejemplo curso fabricii
library(ggthemes)
library(ggalt)
library(maps)
library(rgeos)
library(maptools)
remotes::install_github("hrbrmstr/albersusa")
library(albersusa)
library(grid)
library(sf)
library(terra)
library(dplyr)
library(spData)
library(stars)
library(nngeo)

# Let's load some maps:
#esto crea un vector con los estados de EEUU
states<-usa_composite()  #create a state map thing
smap<-fortify(states,region="fips_state")
counties <- counties_composite()   #create a county map thing
plot(counties)

counties@data <- left_join(counties@data, d2015, by = "fips")
#necesito un vector que tenga cada pixel en un poligono y que cada poligono tenga los valores de pendiente y HFP

slp<-rast("D:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_slope_MX.tif")
hfp<-rast("D:/Maestria_DD/Shapes_MSc_DD/HFP/Dryad_2020/Human_footprint_maps/hfp2013_merisINT_MX.tif")
crp<-vect("D:/Maestria_DD/Shapes_MSc_DD/borrar/extend_prueb.shp")

hfp.res<- hfp |> 
  resample(slp, method= "near") |> 
  writeRaster("./outputs/HFP_MX.tif")

slp.vec <- slp |> 
  as.polygons(dissolve = F) |> 
  st_as_sf() 

hfp.vect <- hfp.res |> 
  as.points() |> 
  st_as_sf()

hfpxslp<-st_join(slp.vec, hfp.vect)

d=expand.grid(x=1:3,y=1:3)
names(d)

leg= ggplot(d, aes(x,y))+
  geom_tile(aes(alpha=x+y,fill=atan(y/x)))+
  scale_fill_viridis_c()+
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  coord_equal()

x=classInt::classIntervals(hfpxslp$wc2.1_30s_slope,4,style = "quantile")
y=classInt::classIntervals(hfpxslp$hfp2013_merisINT_MX,4,style = "quantile")

hfpxslp$x = classInt::findCols(x)
hfpxslp$y = classInt::findCols(y)
hfpxslp$alpha = as.character(hfpxslp$x + hfpxslp$y)
hfpxslp$color = as.character(atan(hfpxslp$y/hfpxslp$x))

map = ggplot()+
  geom_sf(data = hfpxslp,aes(fill=color,alpha=alpha),shape=15, size=11,show.legend = F, color = NA)+
  scale_fill_viridis_d()+
  theme_void()

library(patchwork)
map + leg




# probando glms -----------------------------------------------------------

Z.slp.cats<- read.table("./outputs/Z.slp.cats.txt", sep = "\t", dec=".", header=T)

library(lme4)
library(MASS)
library(moments)
library(tidyverse)
library(easystats)

# generate models
m0.glm <- glm(z ~ 1, family = gaussian, data = Z.slp.cats)
m0.lmer = lmer(z ~ 1 + (1|ecoreg), REML = T, data = Z.slp.cats)

#Testing Random Effects
AIC(logLik(m0.glm))
AIC(logLik(m0.lmer))
#The inclusion of a random effect structure with random intercepts is justified as the AIC of the model with random intercepts is substantially lower than the AIC of the model without random intercepts.


# generate models with 2 different random effect structures
ma.lmer = lmer(z ~ slope*cat + (1|ecoreg), REML = T, data = Z.slp.cats)
mb.lmer = lmer(z ~ slope*cat + (1 + cat | ecoreg), REML = T, data = Z.slp.cats)
mc.lmer = lmer(z ~ slope*cat + (slope*cat||ecoreg), REML = T, data = Z.slp.cats)
md.lmer = lmer(z ~ slope*cat + (1 + slope|ecoreg), REML = T, data = Z.slp.cats)
me.lmer = lmer(z ~ slope*cat + (1 + cat|ecoreg), REML = T, data = Z.slp.cats)
mf.lmer = lmer(z ~ slope*cat + (cat||ecoreg), REML = T, data = Z.slp.cats)

# compare models
anova(ma.lmer, mb.lmer, mc.lmer,md.lmer,me.lmer, test = "Chisq", refit = F)
anova(ma.lmer, mb.lmer, mc.lmer,md.lmer,me.lmer, test = "Chisq", refit = T)




#--------------------------------------------------------------------------------------------------#

#------------------------ hay correlacion entre la pendiente y las bios?! -------------------------

#--------------------------------------------------------------------------------------------------#

spp.dir<-list.files("E:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD", pattern = ".shp$", full.names = T)
spp.list<-list.files("E:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD", pattern = ".shp$", full.names = F) #solo tiene la lista de las especies

#cargar variables
vbles_stack.mx<-rast(c(list.files(path = "./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_MX/", full.names = TRUE)[-c(18,19,10,11)]))

slope.mx<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")

stack<-c(vbles_stack.mx,slope.mx)
names(stack)
list.spp<-list()

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

for (x in 1:length(spp.dir)) {
  #Abrir shapes
  points.sp<-vect(spp.dir[x]) 
  #extraer valores de los rasters
  sp.ext<-extract(x=stack, y=points.sp, ID=FALSE)
  #estandarizar valores extraidos
  sp.ext.sc<-as.data.frame(scale(sp.ext))
  #correlación
  cor_sp<- cor(na.omit(sp.ext.sc), method = "spearman") #sp.ext.sc contiene los valores por pixel de cada punto de presencia estandarizados
  #write matriz
  # write.table(cor_sp, paste("C:/Users/Mauricio Diaz/Documents/Maestria/spp. records_PC/specialists_PC/spec.corr/", substr(spp.list[x],1,nchar(spp.list[x])-4), "_corr.txt", sep=""), sep = "\t", dec=".")
  # #filtro
  cor_sp.long<-flattenCorrMatrix(cor_sp)
  cor_sp.long.filter<-cor_sp.long[which(cor_sp.long[,3]>=-0.8 & cor_sp.long[,3]<=0.8),]
  cor_sp.long.filter$cor<-round(cor_sp.long.filter$cor,3) #con esta función escojo solo los 3 primeros decimales
  cor_sp.long.filter$spp<-spp.list[x]
  list.spp[[x]]<-cor_sp.long.filter
  #write filtro
  # write.table(cor_sp.long.filter, paste("C:/Users/Mauricio Diaz/Documents/Maestria/spp. records_PC/specialists_PC/spec.corr/", substr(spp.list[x],1,nchar(spp.list[x])-4), "_filter.txt", sep=""), sep = "\t", dec=".", row.names=F) 
  #substr(spp.dir[x],1,nchar(spp.dir[x])-4) elimina los últimos 4 caracteres de cada especie (es decir= ".shp")
  
}

list.spp2<-as.data.frame(do.call(rbind, list.spp))
list.spp2 |> 
  as_tibble() |> 
  filter(column=="slope_mx_g",
         cor==-0.323)

