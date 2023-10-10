#Este script tiene como finalidad construir un grafico Range-Diversity a partir de rasters de distribución potencial. Además para colorear los puntos segun una tercera variable

#1. Convertir todos los mapas de distribución en un mapa de presencia-ausencia (PAM)

library(terra) # rast as.data.frame head as.matrix plot summary unique as.factor image t points as.raster rev text diff match extract nrow mean lines
library(sf) # No used functions found
library(tidyverse) # No used functions found # No used functions found # No used functions found
#esto es porque como el objeto CS se demora tanto, lo corri una vez y lo guardé como un objeto .RData
load("F:/Maestria_DD/Shapes_MSc_DD/borrar/PAM_CS_MEX.RData") 

#Direcciones de los rasters de distribución potencial
ras.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = T, pattern =".tif$")
#Lista de nombres de los archivos
ras.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = F, pattern =".tif$") |> 
  substr(1,7) #Esto solo elimina los caracteres .tif

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

#convierto la lista en un DF
spd.list.join<- as.data.frame(do.call(rbind, spd.list), row.names = NULL) 
# write.table(spd.list.join,"F:/Maestria_DD/spp.records_DD/specialists_DD/temporales/spp_PAM.txt", sep="\t", dec=".")
head(spd.list.join)

#construcción de la matriz de presencia ausencia con el paquete letsR
library(letsR) # lets.presab.points
PAM<-lets.presab.points(xy= as.matrix(spd.list.join[,1:2]),
                        species= spd.list.join$spp,
                        xmn = -120,
                        xmx = -85,
                        ymn = 14,
                        ymx = 32,
                        resol = 0.08333)  
x11()
plot(PAM)
summary(PAM)

# ------------------------------------------------------------------------------
# Este script de abajo contiene datos obtenidos de https://github.com/jsoberon/PAMs-Mexico y de un articulo de Sobertón et al. 2021. Contiene además una funcion para realizar el plot de Villalobos et al 2011, el cual contiene dos aproximaciones (sitios y especies) y en este caso son para sitios, en donde se asume que un sitio es un pixel.

# ------------------------------------------------------------------------------
# Packages 
## installing remotes to install biosurvey using it
install.packages("remotes")
library(remotes) # install_github # install_github

## installing biosurey (you may be asked to update other packages; not updating 
## those packages may prevent potential errors)
remotes::install_github("claununez/biosurvey")

## loading other packages
library(biosurvey) # PAM_indices prepare_PAM_CS plot_PAM_CS assign_blocks 
library(maps) # map 

#install.packages("viridis")
library(viridis) # magma viridis viridis_pal 

# function to produce old diagrams -  Plots de Arita et al. 2008, Villalobos et al. 2011
source("https://raw.githubusercontent.com/jsoberon/PAMs-Mexico/master/Function_old_diagrams.R")

# Mi matriz de presencia ausencia -----------------------------------------
PAM.mix<-PAM$Presence_and_Absence_Matrix

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Analysis 
## calculating all indices from the PAM. This is for visualization purposes.
PAM.mix_i <- PAM_indices(PAM = PAM.mix, indices = "all", exclude_column = 1:2)


## performing calculations to produce Christen-Soberon diagrams and their 
## geographic representations. Depending on your computer this process could be 
## time consuming (to make it faster you can use fewer 'randomization_iterations')
PAM.mix_CS <- prepare_PAM_CS(PAM = PAM.mix, exclude_column = 1:2, significance_test = TRUE, 
                             randomization_iterations = 500, CL = 0.05, 
                             keep_randomizations = TRUE, parallel = TRUE)

## the list produced contains all elements needed to produce the following plots

# saving results as RData - Guardar los datos en .RData, porque el objeto CS tarda mucho en construirse
save(PAM.mix, PAM.mix_CS, file = "D:/Mauro/PAM_CS_MEX.RData")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Plotting the PAM and its geographic representation

## colors for richness
### in case some cells were discarded for lack of info, let's match these two
### before defining colors
rich <- PAM.mix_i$Richness_normalized
cols <- magma(length(unique(rich)))
colfact <- as.factor(rich)

## colors for dispersion field
### same matching before defining colors
disf <- PAM.mix_i$Dispersion_field
cols1 <- viridis(length(unique(disf)))
colfact1 <- as.factor(disf)

## plot configuration
x11()
layout(matrix(c(1, 1:3), nrow = 2, byrow = T))

## plotting part of a PAM
par(mar = c(0.5, 2.5, 3, 0.5))
image(t(as.matrix(iucn[, -(1:2)])), axes = FALSE)

axis(2, at = 0.5, labels = "Cells", tick = FALSE)
axis(3, at = 0.5, labels = "Species", tick = FALSE)
axis(3, at = 0.1, labels = "Presence-absence matrix\n", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species richness normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols[colfact], pch = 19)
legend_image <- as.raster(matrix(rev(cols), ncol = 1))
text(x = -115, y = 22, labels = "Richness")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
axis(3, at = -111, labels = "Geographic representation", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species dispersion field normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols1[colfact1], pch = 19)
legend_image <- as.raster(matrix(rev(cols1), ncol = 1))
text(x = -115, y = 22, labels = "Dispersion field")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Christen-Soberon diagrams vs previous diagram
x11()
par(mfrow = c(1, 2))
par(mar = c(4.5, 4.5, 3.5, 0.5))

## original diagram - hice un cambio para poder cambiarle el color
library(grDevices) # x11 as.raster jpeg dev.off rainbow # x11 as.raster jpeg dev.off rainbow
library(tidyverse) # No used functions found # No used functions found # No used functions found
library(viridis) # magma viridis viridis_pal # magma viridis viridis_pal # magma viridis viridis_pal # magma viridis viridis_pal
slope<-rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
cords.pam<-PAM.mix[,1:2]

#1. Extraer valores de pendiente
xy.slope<-terra::extract(slope, PAM.mix[,1:2])
#2. Agregar columna para unir
xy.slope2<- xy.slope |>
  mutate(lon=PAM.mix[,1],
         lat=PAM.mix[,2],
         lonlat=paste0(PAM.mix[,1],PAM.mix[,2])) |> 
  #2.1 Reordenar
  arrange(desc(slope_mx_g)) |> 
  #2.2 columna de color
  mutate(cols=viridis_pal(option = "D")(nrow(xy.slope)))

xy.pam<-PAM.mix[,1:2] |> 
  as.data.frame() |> 
  rename(lon='Longitude(x)',
         lat='Latitude(y)') |> 
  mutate(lonlat=paste0(PAM.mix[,1],PAM.mix[,2]))

xy.merge<- inner_join(xy.pam, xy.slope2, by="lonlat")
head(xy.merge)
cols<-xy.merge$cols

rdp2<-function(mat, view, limits, cols) {
  # First, removes empty columns or rows
  mat=as.matrix(mat)
  j=which(colSums(mat)>0)
  mat2=mat[,j]
  i=which(rowSums(mat2)>0)
  mat3=mat2[i,]
  
  if(view==1){
    # View per sites
    n=dim(mat3)[[1]]
    s=dim(mat3)[[2]]
    alfas<<-rowSums(mat3)
    omegas<<-colSums(mat3)
    alfast<<-alfas/s
    omegast<<-omegas/n
    fist<<-mat3%*%omegast
    fistprom<<-fist/alfas
    
    if(limits==1) {xl=1;yl=1} else {xl=max(fistprom)*1.1;yl=max(alfast)*1.1}
    
    betty<-1/mean(alfast)
    rho<-alfast*(fistprom-1/betty)  
    xM<-seq(1.01/betty,1,length=100)
    yM<-max(rho)/(xM-1/betty)
    xm<-seq(0,.99/betty,length=50)
    ym<-min(rho)/(xM-1/betty)
    
    plot(fistprom,alfast,xlim=c(0,xl),ylim=c(0,yl),
         col = cols,
         xlab="Mean normalized dispersion field",ylab="Normalized richness")
    
    abline(v=1/betty, lty = 2)
    abline(h = min(alfast))
    abline(h = max(alfast))
    lines(xM,yM)
    lines(xM,ym)
  } else {
    #View per species
    mat3t=t(mat3)
    nt=dim(mat3t)[[1]]
    st=dim(mat3t)[[2]]
    alfast<-rowSums(mat3t)
    omegast<-colSums(mat3t)
    alfastt<-alfast/st
    omegastt<-omegast/nt
    fistt<-mat3t%*%omegastt
    fistpromt<-fistt/alfast
    
    # Type of limits: 1 -> limits form 0 to 1, not 1, limits near the maximum of the range of x and y.
    if(limits==1) {xl=1;yl=1} else {xl=max(fistpromt)*1.1;yl=max(alfastt)*1.1}
    
    bettyt<-1/mean(alfastt)
    rhot<-alfastt*(fistpromt-1/bettyt)
    xM<-seq(1.01/bettyt,1,length=100)
    yM<-max(rhot)/(xM-1/bettyt)
    xm<-seq(0,.99/bettyt,length=50)
    ym<-min(rhot)/(xm-1/bettyt)
    
    plot(fistpromt,alfastt,xlim=c(0,xl),ylim=c(0,yl),
         col = cols,
         xlab="Normalized richness",ylab="Mean normalized dispersion field")
    
    abline(v=1/bettyt)
    lines(xM,yM)
    lines(xM,ym)
  }
}

par(cex = 0.6)
rdp2(mat = PAM.mix[, -(1:2)], view = 1, limits = 1, cols=cols)
title(main = "Previous range-diversity plot")

## simple Christen Soberon plot  
par(cex = 0.8)
plot_PAM_CS(PAM.mix_CS, main = "Range-diversity plot", col_all = cols)

#Fin del script