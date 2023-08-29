##### Convertir todos los mapas de distribución en un mapa de presencia-ausencia (PAM)

library(terra)
library(sf)
library(tidyverse)
load("F:/Maestria_DD/Shapes_MSc_DD/borrar/PAM_CS_MEX.RData")
ras.dir<-list.files("D:/Mauro/MB_Mex", full.names = T, pattern =".tif$")
ras.list<-list.files("D:/Mauro/MB_Mex", full.names = F, pattern =".tif$") |> 
  substr(1,7)

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
    #elimino la columna de presencia de los rasters
    dplyr::select(-3)
  names(spd.list)[i]<-ras.list[i]
  print(ras.list[i])
}

spd.list.join<- as.data.frame(do.call(rbind, spd.list), row.names = NULL) 
# write.table(spd.list.join,"F:/Maestria_DD/spp.records_DD/specialists_DD/temporales/spp_PAM.txt", sep="\t", dec=".")

head(spd.list.join)
library(letsR)
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
# Project title: Visualizing Species Richness and Site Similarity from 
#                Presence-absence Matrices
# Authors: Jorge Soberon, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 05/03/2021 (dd/mm/yyyy)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Details:
#
# This script serves to replicate the analyses of for the paper ""Visualizing 
# Species Richness and Site Similarity from Presence-absence Matrices"
# 
# This script and the data used for analysis is stored in a GitHub repository:
# https://github.com/jsoberon/PAMs-Mexico
# 
# The data is a presence-absence matrix (PAM) for 1595 species of terrestrial 
# vertebrates of Mexico. Mexico was divided in a grid of 711 cells of 0.5 decimal
# degrees for purposes of representation. Species ranges used to construct the  
# PAM were obtained from the IUCN database. The number of species remaining 
# after excluding species with no presences in the PAM was 1573.
#
# The analyses presented here are performed using base functions of R, and 
# other functions from the packages biosurvey, maps, and viridis (for colors).
# 
# biosurvey is a new R package available from GitHub (not yet in CRAN), 
# instructions for installing this package are provided below.
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Packages 
## installing remotes to install biosurvey using it
install.packages("remotes")
library(remotes)

## installing biosurey (you may be asked to update other packages; not updating 
## those packages may prevent potential errors)
remotes::install_github("claununez/biosurvey")

## loading other packages
library(biosurvey)
library(maps)

#install.packages("viridis")
library(viridis)

# function to produce old diagrams
source("https://raw.githubusercontent.com/jsoberon/PAMs-Mexico/master/Function_old_diagrams.R")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Data for analysis
## setting directory
setwd("R/Christen-Soberon") # change it to your working directory

## this is the presence absence matrix  
pam_url <- "https://github.com/jsoberon/PAMs-Mexico/blob/master/PAM_IUCN_MEX.RData?raw=true"
load(url(pam_url))

iucn[1:6, 1:8]

# Mi matriz de presencia ausencia -----------------------------------------
PAM.mix<-PAM$Presence_and_Absence_Matrix

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Analysis 
## calculating all indices from the PAM. This is for visualization purposes.
iucn_i <- PAM_indices(PAM = PAM.mix, indices = "all", exclude_column = 1:2)


## performing calculations to produce Christen-Soberon diagrams and their 
## geographic representations. Depending on your computer this process could be 
## time consuming (to make it faster you can use fewer 'randomization_iterations')
iucnCS <- prepare_PAM_CS(PAM = PAM.mix, exclude_column = 1:2, significance_test = TRUE, 
                         randomization_iterations = 500, CL = 0.05, 
                         keep_randomizations = TRUE, parallel = TRUE)

## the list produced contains all elements needed to produce the following plots

# saving results as RData
save(PAM.mix, iucnCS, file = "D:/Mauro/PAM_CS_MEX.RData")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Plotting the PAM and its geographic representation

## colors for richness
### in case some cells were discarded for lack of info, let's match these two
### before defining colors
rich <- iucn_i$Richness_normalized
cols <- magma(length(unique(rich)))
colfact <- as.factor(rich)

## colors for dispersion field
### same matching before defining colors
disf <- iucn_i$Dispersion_field
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

## original diagram
rdp(PAM.mix[, -(1:2)], view = 1, limits = 1)
title(main = "Previous range-diversity plot")

## simple Christen Soberon plot  
plot_PAM_CS(iucnCS, main = "New range-diversity plot")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Christen-Soberon diagram options
colsg <- ifelse(iucnCS$S_significance_id == 0, "gray70", 
                ifelse(iucnCS$S_significance_id == 1, "gray35", "gray1"))

x11()
par(mfrow = c(1, 2), cex = 0.8)

## simple Christen Soberon plot  
plot_PAM_CS(iucnCS, main = "Simple range-diversity plot")

## Christen Soberon plot showing randomized values   
plot_PAM_CS(iucnCS, main = "Random expectations in the plot", col_all = NA,
            add_random_values = T, ylim = c(0, 0.13), add_legend = F)

## Christen Soberon plot showing significant   
plot_PAM_CS(iucnCS, main = "Significant values in the plot", add_significant = T, 
            col_significant_low = "gray35", col_significant_high = "gray1",
            pch_significant_low = 16, pch_significant_high = 16,
            ylim = c(0, 0.13), add_legend = F)

## geographic representation
maps::map(regions = "Mexico")
title(main = "Geographic representation", line = 3.3)
points(PAM.mix[, 1:2], col = colsg, pch = 19)
legend("bottomleft", pch = 19, col = c("gray70", "gray35", "gray1"), 
       legend = c("Non-significant", "Below lower limit", "Above upper limit"), 
       bty = "n", cex = 0.9, inset = -0.02)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Christen-Soberon diagram, exploration in geography
## partitioning diagram area in blocks
n_cols <- 4
n_rows <- 4
data <- data.frame(iucnCS$Richness_normalized, iucnCS$Dispersion_field_normalized)
id <- paste(data[, 1], data[, 2])

xrange <- range(data[, 1])
xinter <- diff(xrange)/n_cols
yrange <- range(data[, 2])
yinter <- diff(yrange)/n_rows
xlb <- seq(xrange[1], xrange[2], xinter)
xlb[length(xlb)] <- xrange[2]
ylb <- seq(yrange[1], yrange[2], yinter)
ylb[length(ylb)] <- yrange[2]
blocks <- assign_blocks(data, 1, 2, n_cols, n_rows, xlb, ylb, block_type = "equal_area")
blocks <- blocks[match(id, paste(blocks[, 1], blocks[, 2])), ]

## sample of 4 blocks to plot
set.seed(1)
sam4 <- sample(unique(blocks$Block), 4)

## block colors 
colsbg1 <- ifelse(blocks$Block == sam4[1], "dodgerblue4", "gray55")
colsbg2 <- ifelse(blocks$Block == sam4[2], "dodgerblue4", "gray55")
colsbg3 <- ifelse(blocks$Block == sam4[3], "dodgerblue4", "gray55")
colsbg4 <- ifelse(blocks$Block == sam4[4], "dodgerblue4", "gray55")


## plotting
x11()
layout(matrix(1:18, nrow = 6), heights = c(1.2, 10, 10, 10, 10, 1), 
       widths = c(1.2, 10, 10))
par(cex = 0.7)

## Christen Soberon diagram in blocks 
### y labels
par(mar = rep(0, 4))
plot.new()
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new()

### title 1
plot.new(); text(0.5, 0.5, labels = "Blocks in range-diversity plot", font = 2, cex = 1)

### CS plots
par(mar = c(2.5, 2, 0.5, 0.5))
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg1)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg2)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg3)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg4)

### x label
par(mar = rep(0, 4))
plot.new(); text(0.5, 0.5, labels = "Normalized richness")

### title 2
plot.new(); text(0.5, 0.5, labels = "Geographic representation", font = 2, cex = 1)

### mapping blocks
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg1, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg2, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg3, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg4, pch = 19)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Exporting figures
## PAM 
jpeg(filename = "Figure_1.jpg", width = 166, height = 130, units = "mm", res = 600)

## plot configuration
layout(matrix(c(1, 1:3), nrow = 2, byrow = T))

## plotting part of a PAM
par(cex = 0.7)
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
text(x = -115, y = 21.5, labels = "Richness")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
axis(3, at = -110.5, labels = "Geographic representation", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species dispersion field normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols1[colfact1], pch = 19)
legend_image <- as.raster(matrix(rev(cols1), ncol = 1))
text(x = -115, y = 21.5, labels = "Dispersion field")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)

dev.off()



## Christen Soberon plots 
jpeg(filename = "./outputs/Figure_3.jpg", width = 166, height = 166, units = "mm", res = 600)

par(mfrow = c(1, 2))
## original diagram - hice un cambio para poder cambiarle el color
library(grDevices)
library(tidyverse)
library(viridis)
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
plot_PAM_CS(iucnCS, main = "Range-diversity plot", col_all = cols)

dev.off()


## Christen Soberon options 
jpeg(filename = "./outputs/Figure_4.jpg", width = 300, height = 166, units = "mm", res = 600)

par(mfrow = c(1, 2), cex = 0.7)

## Christen Soberon plot showing significant   
plot_PAM_CS(iucnCS, main = "Significant values in the plot", add_significant = T, 
            col_significant_low = "gray35", col_significant_high = "gray1",
            pch_significant_low = 15, pch_significant_high = 17,
            ylim = c(0, 0.13), add_legend = T)

## geographic representation
maps::map(regions = "Mexico")
title(main = "Geographic representation", line = 3.3)
points(PAM.mix[, 1:2], col = colsg, pch = 19)
legend("bottomleft", pch = 19, col = c("gray70", "gray35", "gray1"), 
       legend = c("Non-significant", "Below lower limit", "Above upper limit"), 
       bty = "n", cex = 0.9, inset = -0.02)

dev.off()



## Christen Soberon diagram, blocks explorations
jpeg(filename = "Figure_5.jpg", width = 140, height = 200, units = "mm", res = 600)
layout(matrix(1:18, nrow = 6), heights = c(1.2, 10, 10, 10, 10, 1), 
       widths = c(1.2, 10, 10))
par(cex = 0.7)

## Christen Soberon diagram in blocks 
### y labels
par(mar = rep(0, 4))
plot.new()
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new(); text(0.5, 0.57, labels = "Norm. dispersion field / S", srt = 90)
plot.new()

### title 1
plot.new(); text(0.5, 0.5, labels = "Blocks in range-diversity plot", font = 2, cex = 1)

### CS plots
par(mar = c(2.5, 2.5, 0.5, 0.5))
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg1)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg2)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg3)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg4)

### x label
par(mar = rep(0, 4))
plot.new(); text(0.5, 0.5, labels = "Normalized richness")

### title 2
plot.new(); text(0.5, 0.5, labels = "Geographic representation", font = 2, cex = 1)

### mapping blocks
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg1, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg2, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg3, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg4, pch = 19)

dev.off()
# ------------------------------------------------------------------------------





x<-1:20
y<-sample(3:300, size=20)
cols<-rainbow(20)

plot(y,x,col = cols,
     xlab="Mean normalized dispersion field",ylab="Normalized richness",
     pch=4)


a<-data.frame(ID=1:10,let=sample(letters, size=10))
b<-data.frame(ID=sample(1:10, size=10, replace = F),let=sample(LETTERS, size=10))


left_join(a,b)
bind_rows(list(a,b), .id="ID")
inner_join(a,b, by="ID")
b

