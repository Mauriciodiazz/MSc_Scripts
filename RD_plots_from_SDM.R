

#Este script tiene como finalidad construir un grafico Range-Diversity a partir de rasters de distribución potencial. Además para colorear los puntos segun una tercera variable


# Construcción PAM desde raster -------------------------------------------

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
PAM.list.join<- as.data.frame(do.call(rbind, spd.list), row.names = NULL) 
# write.table(PAM.list.join,"F:/Maestria_DD/spp.records_DD/specialists_DD/temporales/spp_PAM.txt", sep="\t", dec=".")
head(PAM.list.join)

PAM.list.join<- read.table("./outputs/tablas/spp_PAM.txt",sep="\t", dec=".", header=T)


# Construccion de PAM -----------------------------------------------------

#construcción de la matriz de presencia ausencia con el paquete letsR
library(letsR) # lets.presab.points
PAM<-lets.presab.points(xy= as.matrix(PAM.list.join[,1:2]),
                        species= PAM.list.join$spp,
                        xmn = -120,
                        xmx = -85,
                        ymn = 14,
                        ymx = 32,
                        resol = 0.08333)  
x11()
plot(PAM)
summary(PAM)


# RD plots  ---------------------------------------------------------------
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
# preparación de datos para el plot de soberon

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


# 1. Plot del RD-plot de Soberon en la geografía -----------------------
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


# RD plots con opcion de cambio de color -------------------------------

# x11()
# par(mfrow = c(1, 2))
# par(mar = c(4.5, 4.5, 3.5, 0.5))


# RD originales Arita y Villalobos ----------------------------------------

# Paquetes 
library(grDevices) # x11 as.raster jpeg dev.off rainbow 
library(tidyverse) # No used functions found 
library(viridis) # magma viridis viridis_pal
library(terra)

# Preparación de voctor de color
slope<-rast("F:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
cat<-rast("./INEGI/Cambios/INEGI_VII_TC.tif")
byclass<-rast("./outputs/temps/by_class_col.tif")
cords.pam<-PAM.mix[,1:2]

slope2<-resample(slope, PAM$Richness_Raster, method="bilinear")
cat2<-resample(cat, PAM$Richness_Raster, method="near")


#1. Extraer valores de pendiente
xy.slope<-terra::extract(c(slope2,cat2,byclass), PAM.mix[,1:2])
head(xy.slope)
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
slp<-xy.merge$slope_mx_g
#clasificacion de color de clasificacion de mapa bivariado
by_class_col<-xy.merge$by_class_col
by_class_col[which(by_class_col==1)]<-"#E7E7E7"
by_class_col[which(by_class_col==2)]<-"#9ed4de"
by_class_col[which(by_class_col==3)]<-"#1bacbe"
by_class_col[which(by_class_col==4)]<-"#f5acac"
by_class_col[which(by_class_col==5)]<-"#b29ea5"
by_class_col[which(by_class_col==6)]<-"#527f8d"
by_class_col[which(by_class_col==7)]<-"#e35959"
by_class_col[which(by_class_col==8)]<-"#ac5255"
by_class_col[which(by_class_col==9)]<-"#5d3f47"
by_class_col[which(is.na(by_class_col))]<-"#000000"


#color de las categorias
cat.col<-xy.merge$INEGI_VII_TC
cat.col[which(cat.col==1)]<-"#248f5d"
cat.col[which(cat.col==2)]<-"#f56038"

# ###### RD plot versión simple -------------------------------------------

## original diagram - hice un cambio para poder cambiarle el color
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
title(main = "Arita & Villalobos range-diversity plot")


# ###### RD plot versión mela ---------------------------------------------

### Function to calculate the RangeDiversity parameters (based on Arita et al. 2012. GEB)

### DATA: 
### Delta = a presence-absence matrix (PAM); ***species are rows and sites are columns*** (without coordinates or any other information)
### *Note that most PAMs nowadays come as transpose; with species as columns and sites as rows. 
### *If so, you'd need to transpose your PAM before using this function. 

RangeDiversity <- function(Delta){
  
  DeltaT <- t(Delta)
  
  #Calculate parameters of range and richness
  
  Species<- nrow(Delta)
  Quadrats<- ncol(Delta)
  One_S<-as.matrix(seq(1, 1, length=Species))
  One_N<-as.matrix(seq(1, 1, length=Quadrats))
  
  SpeciesRichness<- as.matrix(colSums(Delta))
  RangeSize<- as.matrix(rowSums(Delta))
  
  Fill<- sum(RangeSize)
  Beta<- Species*Quadrats/Fill
  
  RangeMean<- mean(RangeSize)
  RangeMin<- min(RangeSize)
  RangeMax<- max(RangeSize)
  D_Volume<- Delta%*%SpeciesRichness
  RangeRichness<- D_Volume/RangeSize
  SpeciesCovariance<-(Delta%*%DeltaT-RangeSize%*%t(RangeSize)/Quadrats)/Quadrats
  RangeVariance<- RangeSize/Quadrats*(1-RangeSize/Quadrats)
  VarianceOfRanges<-t(RangeSize)%*%RangeSize/Species-(t(One_S)%*%RangeSize/Species)^2
  
  RichnessMean<- mean(SpeciesRichness)
  RichnessMin<- min(SpeciesRichness)
  RichnessMax<- max(SpeciesRichness)
  R_Volume<- DeltaT%*%RangeSize
  SiteRange<- R_Volume/SpeciesRichness
  SitesCovariance<- (DeltaT%*%Delta-SpeciesRichness%*%t(SpeciesRichness)/Species)/Species
  RichnessVariance<- SpeciesRichness/Species*(1-SpeciesRichness/Species)
  VarianceOfRichness<-t(SpeciesRichness)%*%SpeciesRichness/Quadrats-(t(One_N)%*%SpeciesRichness/Quadrats)^2
  
  CovarianceSpecies<- RangeSize*(RangeRichness-RichnessMean)/Quadrats/Species
  CovarianceSites<- SpeciesRichness*(SiteRange-RangeMean)/Species/Quadrats
  
  V_species<- VarianceOfRichness/sum(RangeVariance)
  V_sites<- VarianceOfRanges/sum(RichnessVariance)
  U_species<- sum(abs(RangeSize-RangeMean)/Quadrats)/Species
  U_sites<- sum(abs(SpeciesRichness-RichnessMean)/Species)/Quadrats
  
  list(RD.single.values=data.frame(Species, Quadrats, RichnessMean, RangeMean,Fill, Beta, V_species, V_sites),RD.species.values=data.frame(RangeSize,RangeRichness,RangeVariance,CovarianceSpecies),RD.sites.values=data.frame(SpeciesRichness,SiteRange,RichnessVariance,CovarianceSites),SitesCovariance=SitesCovariance,SpeciesCovariance=SpeciesCovariance)  
  
}
sites.RDplot<- function(pamRD, cols) {
  cols<-cols
  
  Species <- pamRD$RD.single.values$Species
  Quadrats <- pamRD$RD.single.values$Quadrats
  RangeSize <- pamRD$RD.species.values$RangeSize
  SpeciesRichness <- pamRD$RD.sites.values$SpeciesRichness
  SiteRange <- pamRD$RD.sites.values$SiteRange
  RangeMin <- min(pamRD$RD.species.values$RangeSize)
  RangeMax <- max(pamRD$RD.species.values$RangeSize)
  RangeMean <- mean(pamRD$RD.species.values$RangeSize)
  
  
  def.par <-par(no.readonly = TRUE)
  xhist<- hist(SiteRange/Quadrats, breaks = seq(0, 1, .05), plot = FALSE)
  yhist<- hist(SpeciesRichness/Species, breaks = seq(0, 1, .05), plot = FALSE)
  topx<-max(xhist$counts)
  topy<-max(yhist$counts)
  nf<- layout(matrix(c(2, 0, 1, 3), 2, 2, T), c(6, 1), c(1, 6), TRUE)
  layout.show(nf)
  
  #RD plot
  par(mar= c(5, 5, 1, 1))
  plot(SiteRange/Quadrats, SpeciesRichness/Species, xlim = c(0,1), ylim = c(0,1), xlab = "Mean proportional per-site range", ylab = "Proportional species richness", pch=21,bg=cols,cex=1.2, col=alpha(cols,1))
  
  
  # #ISOCOVARIANCE LINES
  # x<- seq(RangeMean/Quadrats, 1, length = 100)
  # for(i in c(.01, .05, .1)) {
  #   lines(x, i/(x - RangeMean/Quadrats), lwd = 2, col = "pink")
  # }
  # y<- seq(0, RangeMean/Quadrats, length = 100)
  # for(j in c(-.01, -.05, -.1)) {
  #   lines(y, j/(y- RangeMean/Quadrats), lwd = 2, col = "pink")
  # }
  #
  #
  x<- seq(RangeMin/Quadrats, RangeMean/Quadrats, length = 100)
  lines(x, (RangeMax-RangeMean)/Quadrats/(RangeMax/Quadrats-x), lwd = 2)
  y<- seq(RangeMean/Quadrats, RangeMax/Quadrats, length = 100)
  lines(y, (RangeMean-RangeMin)/Quadrats/(y-RangeMin/Quadrats), lwd = 2)
  segments(RangeMean/Quadrats, 0, RangeMean/Quadrats, 1, lty = 3, lwd = 1.5)

  #Top histogram
  par(mar=c(0, 5, 1, 1))
  barplot(xhist$counts, axes= FALSE, ylim = c(0, topx), space =0)

  #Side histogram
  par(mar=c(5, 0, 1, 1))
  barplot(yhist$counts, axes = FALSE, space = 0, horiz = TRUE)

  par(def.par)
}

RangeDiversity(t(PAM.mix[,-c(1,2)])) |> 
  sites.RDplot(by_class_col)

# RD plot color=categoria ---------------------------------------------

RangeDiversity(t(PAM$Presence_and_Absence_Matrix[,-c(1,2)])) |> 
  sites.RDplot(cat.col)

# Simple Christen Soberon plot   ------------------------------------------
par(cex = 0.8)
plot_PAM_CS(PAM.mix_CS, main = "Range-diversity plot", col_all = cols)


#Fin del script