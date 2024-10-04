# Analisis riqueza vs pendiente usando modelos lineales autoregresivos
# Packages

library(letsR) # lets.presab.points
library(ncf)
library(sf)
library(spatialreg)
library(spdep)
library(terra)
library(tidyverse)
# Analisis riqueza vs pendiente usando modelos lineales autoregresivos

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

#### SAR para 10 puntos
## NOTA: Esto se corre en el pc del lab
#load("Scripts/codigoGLS/SAR_10k.Rdata")

Z.slp.cats<- read.table("D:/Mauro/Msc/Otros/Z.slp.cats.txt", sep = "\t", dec=".", header=T)

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
  
  list.SAR.min<-list()
  list.SAR.max<-list()
  
  # 5. A bucle to do the six models
  for (i in 1:length(scheme)) {
    start <- Sys.time()
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
    
    list.SAR.min[[i]]<-error_min #save each model into a list
    list.SAR.max[[i]]<-error_max #save each model into a list
    print(Sys.time() - start ) #¿What is the timing of each iteration?
  }
  names(list.SAR.min)<-scheme
  names(list.SAR.max)<-scheme
  #6. returning the results in AIC order and save each model list into a list
  return(list(min=list.SAR.min,max=list.SAR.max,results=arrange(results_g, AIC))) 
}

z.slp.sample.Tot<-
  Z.slp.cats |> 
  sample_n(size=10000,replace=F) # replace=F para que no se repita la selección de las filas, es decir, no se sobreestime la pendiente por celdas repetidas


SAR_res_10k.Tot<-
  z.slp.sample.Tot |> 
  SAR_bucle()

SAR_res_10k.Tot[["results"]] #best model: MAX-W
SAR_res_10k.Tot[["max"]]["W"] #best model

z.slp.sample.C<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  sample_n(size=10000,replace=F)

SAR_res_10k.C<-
  z.slp.sample.C |> 
  SAR_bucle()

SAR_res_10k.C[["results"]] #best model: MAX-W
SAR_res_10k.C[["max"]]["W"] #best model

# #Esto genera muestreos aleatorios i veces
# l.SAR<-list()
# for (i in 1:5) {
# 	print(i)
# 	z.slp.sample<-
# 		Z.slp.cats |> 
# 		filter(cat=="C") |> 
# 		sample_n(size=100,replace=F)
# 	
# 	#l.SAR[[i]]<-
# 		z.slp.sample |> 
# 		SAR_bucle()
# 	print(i)
# }


# Correlograms for OLS and SAR models ------------------------------------------------
# best model total area
SAR_res_10k.Tot[["max"]]["W"] #best model
# best model conserved area
SAR_res_10k.C[["max"]]["W"] #best model

# OLS
lm.g<- lm(z ~ slope,
          data = z.slp.sample.Tot)

lm.c<- lm(z ~ slope,
          data = z.slp.sample.C)

hist(residuals(lm.c))
z.slp.sample.C |> 
  ggplot(aes(x=x,y=y, color=z))+
  geom_point(size=0.4)


# 1. Correlograma global ----------------------------------------------------
# OLS
cor.ols.g <- correlog(z.slp.sample.Tot$x, z.slp.sample.Tot$y,
                      z = residuals(lm.g),
                      na.rm = TRUE,
                      increment = 1,
                      resamp = 1)
# SAR
cor.sar.g <- correlog(z.slp.sample.Tot$x, z.slp.sample.Tot$y,
                      z = residuals(SAR_res_10k.Tot[["max"]]["W"]$W),
                      na.rm = TRUE,
                      increment = 1,
                      resamp = 1)

# 2. Correlograma para C -----------------------------------------------------
# OLS
cor.ols.c <- correlog(z.slp.sample.C$x, z.slp.sample.C$y,
                      z = residuals(lm.c),
                      na.rm = TRUE,
                      increment = 1,
                      resamp = 1)
# SAR
cor.sar.c <- correlog(z.slp.sample.C$x, z.slp.sample.C$y,
                      z = residuals(SAR_res_10k.C[["max"]]["W"]$W),
                      na.rm = TRUE,
                      increment = 1,
                      resamp = 1)


par(mfrow=c(2,2))

plot(cor.ols.g, main="National correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.g, main="National correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.c, main="Conserved areas correlation OLS", ylim=c(-3,1), 
     ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.c, main="Conserved areas correlation SAR", ylim=c(-3,1),
     ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

par(mfrow=c(1,1))

# SAR para los biomas ------------------------------------------------------
Z.slp.cats #riqueza + slope + biomas en C y T

# cuántos pixeles hay por cada bioma?
Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME!="Mangroves") |> 
  summarise(count=n(), .by=BIOME_NAME)

# Deserts & Xeric Shrublands (DXS) - 181229 ----
dxs<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Deserts & Xeric Shrublands") |> 
  sample_n(size=10000,replace=F)
head(dxs)

dxs.list<-
  dxs |> 
  SAR_bucle()

dxs.list[["results"]] #best model: MAX-W
dxs.list[["max"]]["W"] #best model
# Correlograma para dxs -----------------------------------------------------

lm.dxs<-
  lm(z ~ slope,
     data = dxs)

# OLS
cor.ols.dxs<- correlog(dxs$x, dxs$y,
                       z = residuals(lm.dxs),
                       na.rm = TRUE,
                       increment = 1,
                       resamp = 1)
# SAR
cor.sar.dxs <- correlog(dxs$x, dxs$y,
                        z = residuals(dxs.list[["max"]]["W"]$W),
                        na.rm = TRUE,
                        increment = 1,
                        resamp = 1)

plot(cor.ols.dxs, main="dxs correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.dxs, main="dxs correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

# Mediterranean Forests, Woodlands & Scrub (MFWS)- 7576 ----
mfws<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Mediterranean Forests, Woodlands & Scrub")

mfws.list<-
  mfws |> 
  SAR_bucle()

mfws.list[["results"]] #best model: MAX-s
mfws.list[["max"]]["S"] #best model
# Correlograma para mfws -----------------------------------------------------

lm.mfws<-
  lm(z ~ slope,
     data = mfws)

# OLS
cor.ols.mfws<- correlog(mfws$x, mfws$y,
                        z = residuals(lm.mfws),
                        na.rm = TRUE,
                        increment = 1,
                        resamp = 1)
# SAR
cor.sar.mfws <- correlog(mfws$x, mfws$y,
                         z = residuals(mfws.list[["max"]]["S"]$S),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)

plot(cor.ols.mfws, main="MFWS correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.mfws, main="MFWS correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")


# Tropical & Subtropical Dry Broadleaf Forests (TSDBF)  95749 ----
tsdbf<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Tropical & Subtropical Dry Broadleaf Forests") |> 
  sample_n(size=10000,replace=F)
head(tsdbf)

tsdbf.list<-
  tsdbf |> 
  SAR_bucle()

tsdbf.list[["results"]] #best model: MAX-W
tsdbf.list[["max"]]["W"] #best model
# Correlograma para tsdbf -----------------------------------------------------

lm.tsdbf<-
  lm(z ~ slope,
     data = tsdbf)

# OLS
cor.ols.tsdbf<- correlog(tsdbf$x, tsdbf$y,
                         z = residuals(lm.tsdbf),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)
# SAR
cor.sar.tsdbf <- correlog(tsdbf$x, tsdbf$y,
                          z = residuals(tsdbf.list[["max"]]["W"]$W),
                          na.rm = TRUE,
                          increment = 1,
                          resamp = 1)

plot(cor.ols.tsdbf, main="tsdbf correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsdbf, main="tsdbf correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")


# Tropical & Subtropical Coniferous Forests (tscf) 271629 ----
tscf<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Tropical & Subtropical Coniferous Forests") |> 
  sample_n(size=10000,replace=F)
head(tscf)

tscf.list<-
  tscf |> 
  SAR_bucle()

tscf.list[["results"]] #best model: MAX-W
tscf.list[["max"]]["W"] #best model
# Correlograma para tscf -----------------------------------------------------

lm.tscf<-
  lm(z ~ slope,
     data = tscf)

# OLS
cor.ols.tscf<- correlog(tscf$x, tscf$y,
                        z = residuals(lm.tscf),
                        na.rm = TRUE,
                        increment = 1,
                        resamp = 1)
# SAR
cor.sar.tscf <- correlog(tscf$x, tscf$y,
                         z = residuals(tscf.list[["max"]]["W"]$W),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)

plot(cor.ols.tscf, main="tscf correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tscf, main="tscf correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")


# Tropical & Subtropical Moist Broadleaf Forests (TSMBF)  80398 ----
tsmbf<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Tropical & Subtropical Dry Broadleaf Forests") |> 
  sample_n(size=10000,replace=F)
head(tsmbf)

tsmbf.list<-
  tsmbf |> 
  SAR_bucle()

tsmbf.list[["results"]] #best model: MAX-W
tsmbf.list[["max"]]["W"] #best model
# Correlograma para tsmbf -----------------------------------------------------

lm.tsmbf<-
  lm(z ~ slope,
     data = tsmbf)

# OLS
cor.ols.tsmbf<- correlog(tsmbf$x, tsmbf$y,
                         z = residuals(lm.tsmbf),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)
# SAR
cor.sar.tsmbf <- correlog(tsmbf$x, tsmbf$y,
                          z = residuals(tsmbf.list[["max"]]["W"]$W),
                          na.rm = TRUE,
                          increment = 1,
                          resamp = 1)

plot(cor.ols.tsmbf, main="tsmbf correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsmbf, main="tsmbf correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")


#  Tropical & Subtropical Grasslands, Savannas & Shrublands (TSGSS)   1908 ----

tsgss<-
  Z.slp.cats |> 
  filter(cat=="C") |> 
  filter(BIOME_NAME=="Tropical & Subtropical Grasslands, Savannas & Shrublands")
head(tsgss)

tsgss.list<-
  tsgss |> 
  SAR_bucle()

tsgss.list[["results"]] #best model: MAX-W
tsgss.list[["max"]]["W"] #best model
# Correlograma para tsgss -----------------------------------------------------

lm.tsgss<-
  lm(z ~ slope,
     data = tsgss)

# OLS
cor.ols.tsgss<- correlog(tsgss$x, tsgss$y,
                         z = residuals(lm.tsgss),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)
# SAR
cor.sar.tsgss <- correlog(tsgss$x, tsgss$y,
                          z = residuals(tsgss.list[["max"]]["W"]$W),
                          na.rm = TRUE,
                          increment = 1,
                          resamp = 1)

plot(cor.ols.tsgss, main="tsgss correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsgss, main="tsgss correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

# Todos en el mismo plot
par(mfrow=c(3,4))

plot(cor.ols.mfws, main="MFWS correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.mfws, main="MFWS correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.dxs, main="DXS correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.dxs, main="DXS correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.tscf, main="TSCF correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tscf, main="TSCF correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.tsdbf, main="TSDBF correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsdbf, main="TSDBF correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.tsgss, main="TSGSS correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsgss, main="TSGSS correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

plot(cor.ols.tsmbf, main="TSMBF correlation OLS", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")
plot(cor.sar.tsmbf, main="TSMBF correlation  SAR", ylim=c(-3,1), ylab="Correlation", xlab="Distance (Clase of distance)")
abline(h=0, col="red")

# Points graph linear models Tot-C -----------------------------------------------
Z.slp.cats

# Best models

# Total
SAR_res_10k.Tot[["results"]] #max - W
summary(SAR_res_10k.Tot[["max"]]["W"]$W, Nagelkerke=TRUE)

# Conserved
SAR_res_10k.C[["results"]] #max - W
summary(SAR_res_10k.C[["max"]]["W"]$W, Nagelkerke=TRUE)

#summary(error_min_t, Nagelkerke=TRUE) Este ya no va

confint(dxs.list[["max"]]["W"]$W, level=0.95)

slp_z_g<-
z.slp.sample.Tot |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=2.536156, slope=0.471712), color="red") +
  annotate(geom="text", x=20, y=30, 
           label=expression(atop("y=0.47x + 2.54",
                                 paste(italic("p"),"-value<2.2e-16"))))+
  labs(x="Terrain slope", y="Richness", title="National") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

slp_z_cat<-
z.slp.sample.C |> 
  ggplot(aes(y=z, x=slope, color=cat))+
  geom_point()+
  geom_abline(aes(intercept=2.5060748, slope=0.2571645), color="red") +
  annotate(geom="text", x=25, y=28, 
           label=expression(atop("y=0.26x + 2.51",
                           paste(italic("p"),"-value<2.2e-16"))))+
  scale_color_manual(values=c("#248f5d"), labels = c("Conservado"))+
  labs(x="Terrain slope", y="Richness", color="Categorias", title="Conserved") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = "none")

#library(patchwork)
slp_z_g + slp_z_cat + plot_annotation(tag_levels = 'A')

ggsave(filename = "./outputs/SAR_zxslp.png",
       width = 20,
       height = 10, #alto
       scale=1,
       units ="cm",
       dpi = 200)

# Points graph linear models biomes -----------------------------------------------

# DXS: Desiertos y matorrales xerófilos
dxs.list[["results"]] #max-W
dxs.list[["max"]][["W"]] |> summary(Nagelkerke=TRUE)
confint(dxs.list[["max"]]["W"]$W, level=0.95)

plot.dxs<-
  dxs |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=1.4589790, slope=0.1821844), color="red") +
  annotate(geom="text", x=21, y=14.5, 
           label=expression(atop("y=0.18x + 1.46",
                                 paste(italic("p"),
                                       "-value<2.2e-16"))), 
           parse=T, size=3.8)+
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="DXS") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())


# MFWS: Mediterranean forests, woodlands, and scrub
mfws.list[["results"]] #max-S
mfws.list[["max"]][["S"]] |> summary(Nagelkerke = T)
confint(mfws.list[["max"]]["S"]$S, level=0.95)

plot.mfws<-
mfws |> 
  ggplot(aes(y=as.character(z), x=slope))+
  geom_point()+
  geom_abline(aes(intercept=1.01114942, slope=-0.00074703), color="red") +
  annotate(geom="text", x=18, y=2.3, 
           label=expression(atop("y=-0.0007x + 1.01",
                                 paste(italic("p"),"-value=0.34"))), 
           parse=T, size=3.8)+
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="MFWS") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# TSDBF: Tropical and subtropical dry broadleaf forests
tsdbf.list[["results"]] #max-W
tsdbf.list[["max"]][["W"]] |> summary(Nagelkerke = T)
confint(tsdbf.list[["max"]]["W"]$W, level=0.95)

plot.tsdbf<-
  tsdbf |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=5.389625, slope=0.168228), color="red") +
  annotate(geom="text", x=25, y=24, 
           label=expression(atop("y=0.17x + 5.39",
                                 paste(italic("p"),"-value<2.2e-16"))), 
           parse=T, size=3.8)+
  xlim(0, 33) +  ylim(0, 28) +
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="TSDBF") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# TSCF: Bosques de coníferas tropicales y subtropicales
tscf.list[["results"]] #max-W
tscf.list[["max"]][["W"]] |> summary(Nagelkerke = T)
confint(tscf.list[["max"]]["W"]$W, level=0.95)

plot.tscf<-
tscf |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=4.2231219, slope=0.0948695), color="red") +
  annotate(geom="text", x=24, y=25, 
           label=expression(atop("y=0.09x + 4.22",
                                 paste(italic("p"),"-value<2.2e-16"))), 
           parse=T, size=3.8)+
  xlim(0, 30) +  ylim(0, 28) +
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="TSCF") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

# TSMBF: Tropical and subtropical moist broadleaf forests
tsmbf.list[["results"]] #max-W
tsmbf.list[["max"]][["W"]] |> summary(Nagelkerke = T)
confint(tsmbf.list[["max"]]["W"]$W, level=0.95)


plot.tsmbf<-
  tsmbf |> 
  ggplot(aes(y=z, x=slope))+
  geom_point()+
  geom_abline(aes(intercept=4.8305362, slope=0.1640939), color="red") +
  annotate(geom="text", x=25, y=24, 
           label=expression(atop("y=0.16x + 4.83",
                                 paste(italic("p"),"-value<2.2e-16"))),
           parse=T, size=3.8)+
  xlim(0, 32) +  ylim(0, 28) +
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="TSMBF") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

# TSGSS: Tropical and subtropical grasslands
tsgss.list[["results"]] #max-W
tsgss.list[["max"]][["W"]] |> summary(Nagelkerke = T)
confint(tsgss.list[["max"]]["W"]$W, level=0.95)

plot.tsgss<-
tsgss |> 
  ggplot(aes(y=as.character(z), x=slope))+
  geom_point()+
  geom_abline(aes(intercept=2.133589, slope=0.098550), color="red") +
  annotate(geom="text", x=1.5, y=3, 
           label=expression(atop("y=0.1x + 2.13",
                                 paste(italic("p"),"-value=0.01"))),
           parse=T, size=3.8)+
  labs(x="Terrain slope", y="Richness", color="Categorias", 
       title="TSGSS") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

plot.dxs + plot.mfws + plot.tsdbf + plot.tscf + plot.tsmbf + plot.tsgss

ggsave(filename = "./outputs/SAR_biomes.png",
       width = 20,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)

(slp_z_g + slp_z_cat ) / (plot.dxs + plot.mfws + plot.tsdbf + plot.tscf + plot.tsmbf + plot.tsgss) + plot_annotation(tag_levels = list(c('A', ' ', 'B')))

ggsave(filename = "./outputs/SAR_merged.png",
       width = 13,
       height = 15, #alto
       scale=2,
       units ="cm",
       dpi = 200)

#save.image("./Scripts/codigoGLS/SAR_10k.Rdata")

# FIN DEL SCRIPT


