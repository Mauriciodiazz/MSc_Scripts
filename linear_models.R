# Analisis riqueza vs pendiente usando modelos lineales autoregresivos

# Packages
library(letsR) # lets.presab.points
library(ncf)
library(sf)
library(spatialreg)
library(spdep)
library(tidyverse)
library(terra)

# Mapa de riqueza a partir de raster con LetsR -----------------------------------------------
# Esta parte del script realiza una PAM desde rasters, esto con el fin de establecer de manera más sencilla mapas de riqueza de diversa resolución  

#Direcciones de los rasters de distribución potencial
ras.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = T, pattern =".tif$")
#Lista de nombres de los archivos
ras.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = F, pattern =".tif$") |> 
  substr(1,7) #Esto solo elimina los caracteres .tif

# Crear PAM y raster de riqueza con Lets R a partir de rasters ------------

# En este bucle lo que hago es abrir los rasters, los reescalo si es necesario y luego los convierto en puntos, de manera que el centroide del pixel contiene valores de 0 y 1. Esto lo guardo en una lista y cada objeto de la lista contiene la lista de presencias (valores de 1 de los rasters de distribución potencial)

spd.list<-list()
for(i in 1:length(ras.list)){
  spd.list[[i]]<-rast(ras.dir[i]) |> 
    #hay que reescalar el raster? a 30km = 0.4166667
    #  aggregate(fact=10, fun= 'modal') |> 
    #convierto en un df
    as.data.frame(xy=T) |> 
    #le pongo una columna con el nombre de la spp
    mutate(spp=ras.list[i]) |> 
    #selecciono solo los valores de "presencia"
    filter(.data[[names(rast(ras.dir[i]))]]==1) |> 
    #elimino la columna con valores de 1
    dplyr::select(-3)
  names(spd.list)[i]<-ras.list[i]
  print(ras.list[i])
}

## convierto la lista en un DF y lo guardo para optimizar tiempo
spd.list.join<- as.data.frame(do.call(rbind, spd.list), row.names = NULL) 
# write.table(spd.list.join,"F:/Maestria_DD/spp.records_DD/specialists_DD/temporales/spp_PAM.txt", sep="\t", dec=".", row.names=F)
head(spd.list.join)

#Aqui voy a hacer lo mismo de arriba, pero la diferencia es que a cada área de distribución la voy a cortar para las zonas conservadas
cons<-vect("INEGI/Cambios/TC_s7/C_disuelto_s7.shp")

spd.list<-list()
for(i in 1:length(ras.list)){
  spd.list[[i]]<-rast(ras.dir[i]) |> 
    #corto el raster al polígono de áreas conservadas
    crop(cons,mask=T, touches=F, 
         filename=paste0("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/MB_C/",ras.list[i],"_C.tif"),
         overwrite=T) |> 
    #convierto en un df
    as.data.frame(xy=T) |> 
    #le pongo una columna con el nombre de la spp
    mutate(spp=ras.list[i]) |> 
    #selecciono solo los valores de "presencia"
    filter(.data[[names(rast(ras.dir[i]))]]==1) |> 
    #elimino la columna con valores de 1
    dplyr::select(-3)
  names(spd.list)[i]<-ras.list[i]
  print(ras.list[i])
}

## convierto la lista en un DF y lo guardo para optimizar tiempo
as.data.frame(do.call(rbind, spd.list), row.names = NULL) |> 
  write.table("./outputs/tablas/spp_PAM_C.txt", sep="\t", dec=".", row.names=F)
head(spd.list.join)


####

library(tidyterra)
mx.sp<-vect("./Mexico_Estados/Mexico_continent.shp")
spd.list.join[1:2000,] |> 
  ggplot(aes(x=x,y=y)) + 
  geom_point()+
  geom_spatvector(data=mx.sp, fill=NA)

spd.list.join |> 
  filter(x<(-115)) |> 
  nrow()

ggplot(mx.sp) + 
  geom_spatvector()+
  geom_point(data=spd.list.join |> 
               filter(x<(-115)), 
             aes(x=x,y=y))

####
# Cargo la tabla para optimizar tiempo
#spd.list.join contiene las coordenadas de los centroides de los pixeles de distribución de cada especie
spd.list.join<- read.table("./outputs/tablas/spp_PAM.txt",sep="\t", dec=".", header=T)
head(spd.list.join)

#Para poder hacer los modelos de regresión, es necesario crear un mapa de riqueza y con letsR puedo crear un mapa de riqueza apartir de una tabla de las coordanas de las especies (spd.list.join) 

#1. Hay que crear la Matriz de presencia-ausencia PAM a partir de los puntos de presencia a una resolucion de 10km
PAM<-lets.presab.points(xy= as.matrix(spd.list.join[,1:2]),
                        species= spd.list.join$spp,
                        xmn = -120,
                        xmx = -85,
                        ymn = 14,
                        ymx = 33, 
                        resol = 0.5) #10km 0.08333 / 30km 0.24999 /50km 0.41665

x11()
plot(PAM)
summary(PAM)
PAM.raster<-PAM$Richness_Raster

# 2. Cargar y remuestrear el rasters a la resolución de la capa de riqueza ----------

# Remuestrear pendiente ---------------------------------------------------

slope<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
slope.3<-resample(slope, PAM.raster, method="bilinear")

# Remuetsrear áreas conservadas y transformadas --------------

inegivii<-rast("./INEGI/Cambios/INEGI_VII_TC.tif")
inegivii.3<-resample(inegivii, PAM.raster, method="mode")

# 3. Extraer pendiente e inegi vii para cada pix de riqueza --------

z.slp<-terra::extract(c(slope.3,inegivii.3, PAM.raster), 
               as.points(PAM.raster), xy=T) |> 
  na.omit()|> 
  filter(lyr.1!=0)

# ordenar columnas
z.slp<-z.slp[,c(5,6,4,2,3)]
names(z.slp)<-c("x","y","z","slope","cat")
z.slp$cat<-as.character(z.slp$cat)
z.slp$cat2<-z.slp$cat
z.slp[which(z.slp$cat2==2),6]<-"T"
z.slp[which(z.slp$cat2==1),6]<-"C"

str(z.slp)
levels(factor(z.slp$cat2))

#write.table(z.slp, "./outputs/tablas/z.slp_05.txt", sep="\t", dec=".", row.names=F)
z.slp<-read.table("./outputs/tablas/z.slp_05.txt", sep="\t", dec=".", header=T)
head(z.slp)

# OLS ---------------------------------------------------------------------

# Primero, un modelo simple para observar si existe autocorrelación espacial

lm.simple<- lm(z ~ slope*cat2,
               data = z.slp)

cor.lm.res <- correlog(z.slp$x, z.slp$y,
                       z = residuals(lm.simple),
                       na.rm = TRUE,
                       increment = 1,
                       resamp = 1)
plot(cor.lm.res)

### Existe autocorrelación espacial!


# Modelos SAR  ------------------------------------------------------------
# Estos modelos  incorporan la autocorrelación espacial a través de una lista de pesos según la cantidad de vecinos que tenga un pixel. De manera, que asume que aquellos pixeles con muchos vecinos, entonces tendrán un mayor valor de peso y por ende en esos lugares la relación y~x va a ser distinta.

#Este método utiliza 3 tipos de pesos espaciales: c("W", "C", "S") y dos distancias (min y max), por lo tanto, hay que establecer 6 modelos Wmax, Wmin, Cmax, Cmin, etc. Y elegir cuál de ellos es el mejor a través de una evaluación con el AIC.

# preparación de una función para correr con el modelo global y el de cada categoría. Esta función coje la base de datos y realiza los 6 modelos y extrae información importante para cada uno y realiza una tabla

SAR_bucle<-function(data){
  #1. DF to contain the model results and AIC values
  results_g <- 
    data.frame(
      type=rep(c("min","max"), each = 3),
      style=NA,
      NK_Rsquared = NA,
      lamnda = NA,
      LR_ratio=NA,
      LR_pval=NA,
      Ase=NA, #Asymptotic standard error
      Ase_pval=NA,
      Ase_z=NA,
      Wald=NA,
      Wald_pval=NA,
      log_ll=NA,
      sigma_sq=NA,
      AIC = NA,
      AIC_lm= NA,
      Imoran=NA)
  
  # 2. Data spacialization
  sp <- st_as_sf( x = data, 
                  coords = c("x","y"),
                  crs = 4326)
  
  # 3. Create neighbours list of class nb --------------------------------------
  k1 = knn2nb(knearneigh(sp, k = 1))#k hace referencia a los pixeles cercanos de cada pixel. 
  #1 solo cuenta la distancia al de la derecha; 2 derecha-izquierda; 4 costados; 8 todos los cercanos
  
  # 4. Calculate distances -----------------------------------------------------
  dist <- unlist(nbdists(k1, sp, longlat=T)) #cuidado! si sp es una matriz, entonces debe especificarse longlat=T, porque la medida de distancia debe ser en grados
  max1 <- max(dist)
  min1 <- min(dist)
  
  # 4.1 create a series of neighbour matrices based on different distances (distances are in km)
  d_min <- dnearneigh(sp, longlat = F, d1 = 0, d2 = min1) #d1 y d2 debe ser en grados
  d_max <- dnearneigh(sp, longlat = F, d1 = 0, d2 = max1)
  
  # Spatial weigths creation ------------------------------------------------
  scheme <- c("W", "C", "S")
  
  # 5. A bucle to do the six models
  for (i in 1:length(scheme)) {
    # Spatial weights
    sw_min <- nb2listw(d_min, zero.policy = TRUE, style = scheme[i]) 
    sw_max <- nb2listw(d_max, zero.policy = TRUE, style = scheme[i]) 
    
    # minimum and maximum SAR
    error_min <- spatialreg::errorsarlm(
      z ~ slope,
      data = data,
      listw = sw_min,
      tol = 1e-12, 
      zero.policy = TRUE)
    
    error_max <- spatialreg::errorsarlm(
      z ~ slope,
      data = data,
      listw = sw_max,
      tol = 1e-12, 
      zero.policy = TRUE)
    
    # Model's summary
    error_min_s <-summary(error_min, Nagelkerke=TRUE)
    error_max_s <-summary(error_max, Nagelkerke=TRUE)
    
    results_g$style[i]<-scheme[i]
    results_g$style[i+3]<-scheme[i]
    
    results_g$NK_Rsquared[i]<-error_min_s$NK
    results_g$NK_Rsquared[i+3]<-error_max_s$NK
    
    results_g$lamnda[i]<-error_min_s$lambda[[1]]
    results_g$lamnda[i+3]<-error_max_s$lambda[[1]]
    
    results_g$LR_ratio[i]<-error_min_s$LR1$statistic[1]
    results_g$LR_ratio[i+3]<-error_max_s$LR1$statistic[1]
    
    results_g$LR_pval[i]<-error_min_s$LR1$p.value
    results_g$LR_pval[i+3]<-error_max_s$LR1$p.value
    
    results_g$Wald[i]<-error_min_s$Wald1$statistic
    results_g$Wald[i+3]<-error_max_s$Wald1$statistic
    
    results_g$Wald_pval[i]<-error_min_s$Wald1$p.value
    results_g$Wald_pval[i+3]<-error_max_s$Wald1$p.value
    
    results_g$log_ll[i]<-error_min_s$LL
    results_g$log_ll[i+3]<-error_max_s$LL
    
    results_g$sigma_sq[i]<-error_min_s$s2
    results_g$sigma_sq[i+3]<-error_max_s$s2
    
    results_g$AIC[i]<- AIC(error_min)
    results_g$AIC[i+3]<- AIC(error_max)
    
    results_g$AIC_lm[i]<- error_min_s$AIC_lm.model
    results_g$AIC_lm[i+3]<- error_max_s$AIC_lm.model
    
    results_g$Imoran[i]<- moran(data$z,  sw_min, length(k1), Szero(sw_min))$I
    results_g$Imoran[i+3]<- moran(data$z,  sw_max, length(k1), Szero(sw_max))$I
    
  }
  
  #6. returning the results in AIC order
  return(results_g |> 
           arrange(AIC)) 
  
}

# Como ya se cuales son los mejores modelos, los voy a correr cada uno para los graficos y las tablas

# Modelo Global -----------------------------------------------------------

lm.g<- lm(z ~ slope,
               data = z.slp)

SAR_bucle(z.slp)
#Mejor modelo global= max_w


sp_g <- st_as_sf( x = z.slp, 
                coords = c("x","y"),
                crs = 4326)

k1_g = knn2nb(knearneigh(sp_g, k = 1))


dist_g <- unlist(nbdists(k1_g, sp_g, longlat = T))

max1_g <- max(dist_g)

d_max_g <- dnearneigh(sp_g, longlat = F, d1 = 0, d2 = max1_g)

sw_max_g <- nb2listw(d_max_g, zero.policy = TRUE, style = "W") 

# SAR
error_max_g <- spatialreg::errorsarlm(
  z ~ slope,
  data = z.slp,
  listw = sw_max_g,
  tol = 1e-12, 
  zero.policy = TRUE)

summary(error_max_g, Nagelkerke=TRUE)

plot(z.slp$slope, z.slp$z)

# Modelo para C -----------------------------------------------------------
C<-
  z.slp |> 
  filter(cat2=="C")

lm.c<- lm(z ~ slope,
               data = C)

SAR_bucle(C)

#Mejor modelo para C= min_w

sp_c <- st_as_sf( x = C, 
                  coords = c("x","y"),
                  crs = 4326)

k1_c = knn2nb(knearneigh(sp_c, k = 1)) 


dist_c <- unlist(nbdists(k1_c, sp_c, longlat = T))

max1_c <- max(dist_c)

d_max_c <- dnearneigh(sp_c, longlat = F, d1 = 0, d2 = max1_c)

sw_max_c <- nb2listw(d_max_c, zero.policy = TRUE, style = "W") 

# SAR
error_max_c <- spatialreg::errorsarlm(
  z ~ slope,
  data = C,
  listw = sw_max_c,
  tol = 1e-12, 
  zero.policy = TRUE)

summary(error_max_c, Nagelkerke=TRUE)

z~slope*cat2

# 
# # Modelo para T -----------------------------------------------------------
# Tr<-
#   z.slp |> 
#   filter(cat2=="T")
# 
# lm.t<- lm(z ~ slope,
#           data = Tr)
# 
# SAR_bucle(Tr)
# 
# #Mejor modelo para C= max_w
# 
# sp_t <- st_as_sf( x = Tr, 
#                   coords = c("x","y"),
#                   crs = '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')
# 
# k1_t = knn2nb(knearneigh(sp_t, k = 1)) 
# 
# 
# dist_t <- unlist(nbdists(k1_t, sp_t))
# 
# min1_t <- min(dist_t)
# 
# d_min_t <- dnearneigh(sp_t, longlat = F, d1 = 0, d2 = min1_t)
# 
# sw_min_t <- nb2listw(d_min_t, zero.policy = TRUE, style = "W") 
# 
# # SAR
# error_min_t <- spatialreg::errorsarlm(
#   z ~ slope,
#   data = Tr,
#   listw = sw_min_t,
#   tol = 1e-12, 
#   zero.policy = TRUE)
# 
# summary(error_min_t, Nagelkerke=TRUE)


# Correlograms for OLS and SAR models ------------------------------------------------
error_max_g
error_max_c
#error_min_t ya no va 

lm.g
lm.c
#lm.t ya no va 


# correlograma global ----------------------------------------------------
# OLS
cor.ols.g <- correlog(z.slp$x, z.slp$y,
                         z = residuals(lm.g),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)
# SAR
cor.sar.g <- correlog(z.slp$x, z.slp$y,
                         z = residuals(error_min_g),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)

# Correlograma para C -----------------------------------------------------
# OLS
cor.ols.c <- correlog(C$x, C$y,
                         z = residuals(lm.c),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)
# SAR
cor.sar.c <- correlog(C$x, C$y,
                         z = residuals(error_max_c),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)

# # Correlograma para T ---------------------------------------------------
# # OLS
# cor.ols.t <- correlog(Tr$x, Tr$y,
#                          z = residuals(lm.t),
#                          na.rm = TRUE,
#                          increment = 1,
#                          resamp = 1)
# # SAR
# cor.sar.t <- correlog(Tr$x, Tr$y,
#                          z = residuals(error_min_t),
#                          na.rm = TRUE,
#                          increment = 1,
#                          resamp = 1)

par(mfrow=c(2,2))

plot(cor.ols.g, main="Correlación global OLS", ylim=c(-3,1), ylab="Correlación", xlab="Distancia (Clase de distancia)")
abline(h=0, col="red")
plot(cor.sar.g, main="Correlación global SAR", ylim=c(-3,1), ylab="Correlación", xlab="Distancia (Clase de distancia)")
abline(h=0, col="red")

plot(cor.ols.c, main="Correlación áreas conservadas OLS", ylim=c(-3,1), 
     ylab="Correlación", xlab="Distancia (Clase de distancia)")
abline(h=0, col="red")
plot(cor.sar.c, main="Correlación áreas conservadas SAR", ylim=c(-3,1),
     ylab="Correlación", xlab="Distancia (Clase de distancia)")
abline(h=0, col="red")

# plot(cor.ols.t, main="OLS transformed areas correlation", ylim=c(-3,1))
# abline(h=0, col="red")
# plot(cor.sar.t, main="SAR transformed areas correlation", ylim=c(-3,1))
# abline(h=0, col="red")

par(mfrow=c(1,1))



# Points graph linear models -----------------------------------------------

summary(error_max_g, Nagelkerke=TRUE)
summary(error_max_c, Nagelkerke=TRUE)
#summary(error_min_t, Nagelkerke=TRUE) Este ya no va

slp_z_g<-
z.slp |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=6.422968, slope=1.656030), color="#000000") +
  annotate(geom="text", x=2, y=37, label="y = 1.66x + 6.42", color="black", fontface = 'bold')+
  labs(x="Pendiente del terreno", y="Riqueza de especies") +
  theme_classic()

slp_z_cat<-
z.slp |> 
  filter(cat2=="C") |> 
  ggplot(aes(y=z, x=slope, color=cat2))+
  geom_point()+
  geom_abline(aes(intercept=4.70728, slope=1.87103), color="#248f5d") +
#  geom_abline(aes(intercept=3.85030, slope=1.68200), color="#e5b636") +
  annotate(geom="text", x=2, y=37, label="y = 1.87x + 4.70", color="#248f5d", fontface = 'bold')+
 # annotate(geom="text", x=3, y=35, label="y = 1.68 x + 3.85", color="#e5b636")+
  scale_color_manual(values=c("#248f5d"), labels = c("Conservado"))+
  labs(x="Pendiente del terreno", y="Riqueza de especies", color="Categorias") +
  theme_classic()+
  theme(axis.title.y = element_blank(),
        legend.position = "none")

#library(patchwork)
slp_z_g + slp_z_cat + plot_annotation(tag_levels = 'A')

ggsave(filename = "./outputs/SAR_zxslp.svg",
       width = 20,
       height = 10, #alto
       scale=1,
       units ="cm",
       dpi = 200)

# FIN DEL SCRIPT


