# Mapa bivariado

library(raster)
library(readxl)
library(stringr)
library(bivariatemaps)
library(dplyr)

#elev<-raster("D:/AYUDANTIA_LAB_BIOCLIMATOLOGÍA/Elevacion_slope/elevation_1KMmn_GMTEDmn.tif")
slp<-raster("D:/Maestria_DD/Shapes_MSc_DD/WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_slope_MX.tif")
#hft<-raster("D:/Proyectos/En proceso/Pendiente_deforestacion/Bivariados/HFP_MX.tif")
hft<-raster("D:/Maestria_DD/Shapes_MSc_DD/HFP/Dryad_2020/Human_footprint_maps/hfp2013_merisINT_MX.tif")

#esto es para calcular la pendiente
slope2<-crop(slope, hft, snap="near")
hft
slope3<-resample(slope2, hft) %>% mask(hft)
slp<- terrain(elev3, v="slope", unit="degrees")
plot(hft)

vect("C:/Users/Mauricio Diaz/Downloads/borrar/MEX_ADM0.shp/MEX_ADM0.shp") |> 
  plot()

mex<- sf::read_sf("D:/Maestria_DD/Shapes_MSc_DD/Mexico_Estados/Mexico_continent.shp")

hft.mx<-mask(hft, mex)
slp.mx<-mask(slp, mex)

hft.p<- rasterToPoints(hft.mx) %>% as.data.frame()
summary(hft.p)
hft.p<- hft.p %>% filter(hfp2013_merisINT_MX<100)

hft.slp<-extract(slp.mx, hft.p[,1:2]) %>% cbind(hft.p)
head(hft.slp)

names(hft.slp)[1]<-"Slope"
hft.slp2<-na.omit(hft.slp)

slope2<-rasterize(hft.slp2[,2:3], slp.mx, field=hft.slp2$Slope, silent=T)
plot(slope2)

hfp2<-rasterize(hft.slp2[,2:3], slp.mx, field=hft.slp2$hfp2013_merisINT_MX)
plot(hfp2)

# Define the number of breaks
nBreaks <- 8

# Create the colour matrix first the function

colmat <- function(nbreaks = 3, breakstyle = "quantile",
                   upperleft = "#0096EB", upperright = "#820050", 
                   bottomleft = "#BEBEBE", bottomright = "#FFE60F",
                   xlab = "x label", ylab = "y label", plotLeg = TRUE,
                   saveLeg = FALSE) {
    # TODO - replace any tidyr, dplyr etc. functions with data.table #
    library(tidyverse)
    require(ggplot2)
    require(classInt)
    if (breakstyle == "sd") {
        warning("SD breaks style cannot be used.\nWill not always return the correct number of breaks.\nSee classInt::classIntervals() for details.\nResetting to quantile",
                call. = FALSE, immediate. = FALSE)
        breakstyle <- "quantile"}
    # The colours can be changed by changing the HEX codes for:
    # upperleft, upperright, bottomleft, bottomright
    # From http://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/
    # upperleft = "#64ACBE", upperright = "#574249", bottomleft = "#E8E8E8", bottomright = "#C85A5A",
    # upperleft = "#BE64AC", upperright = "#3B4994", bottomleft = "#E8E8E8", bottomright = "#5AC8C8",
    # upperleft = "#73AE80", upperright = "#2A5A5B", bottomleft = "#E8E8E8", bottomright = "#6C83B5", 
    # upperleft = "#9972AF", upperright = "#804D36", bottomleft = "#E8E8E8", bottomright = "#C8B35A",
    # upperleft = "#DA8DC8", upperright = "#697AA2", bottomleft = "#E8E8E8", bottomright = "#73BCA0",
    # Similar to Teuling, Stockli, Seneviratnea (2011) [https://doi.org/10.1002/joc.2153]
    # upperleft = "#F7900A", upperright = "#993A65", bottomleft = "#44B360", bottomright = "#3A88B5",
    # Viridis style
    # upperleft = "#FEF287", upperright = "#21908D", bottomleft = "#E8F4F3", bottomright = "#9874A1",
    # Similar to Fjeldsa, Bowie, Rahbek 2012
    # upperleft = "#34C21B", upperright = "#FFFFFF", bottomleft = "#595757",  bottomright = "#A874B8",
    # Default from original source
    # upperleft = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", bottomright = "#FFE60F",
    my.data <- seq(0, 1, .01)
    # Default uses terciles (Lucchesi and Wikle [2017] doi: 10.1002/sta4.150)
    my.class <- classInt::classIntervals(my.data,
                                         n = nbreaks,
                                         style = breakstyle,
    )
    my.pal.1 <- classInt::findColours(my.class, c(upperleft, bottomleft))
    my.pal.2 <- classInt::findColours(my.class, c(upperright, bottomright))
    col.matrix <- matrix(nrow = 101, ncol = 101, NA)
    for (i in 1:101) {
        my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
        col.matrix[102 - i, ] <- classInt::findColours(my.class, my.col)
    }
    ## need to convert this to data.table at some stage.
    col.matrix.plot <- col.matrix %>%
        as.data.frame(.) %>% 
        mutate("Y" = row_number()) %>%
        mutate_at(.tbl = ., .vars = vars(starts_with("V")), .funs = list(as.character)) %>% 
        pivot_longer(data = ., cols = -Y, names_to = "X", values_to = "HEXCode") %>% 
        mutate("X" = as.integer(sub("V", "", .$X))) %>%
        distinct(as.factor(HEXCode), .keep_all = TRUE) %>%
        mutate(Y = rev(.$Y)) %>% 
        dplyr::select(-c(4)) %>%
        mutate("Y" = rep(seq(from = 1, to = nbreaks, by = 1), each = nbreaks),
               "X" = rep(seq(from = 1, to = nbreaks, by = 1), times = nbreaks)) %>%
        mutate("UID" = row_number())
    # Use plotLeg if you want a preview of the legend
    if (plotLeg) {
        p <- ggplot(col.matrix.plot, aes(X, Y, fill = HEXCode)) +
            geom_tile() +
            scale_fill_identity() +
            coord_equal(expand = FALSE) +
            theme_void() +
            theme(aspect.ratio = 1,
                  axis.title = element_text(size = 12, colour = "black",hjust = 0.5, 
                                            vjust = 1),
                  axis.title.y = element_text(angle = 90, hjust = 0.5)) +
            xlab(bquote(.(xlab) ~  symbol("\256"))) +
            ylab(bquote(.(ylab) ~  symbol("\256")))
        print(p)
        assign(
            x = "BivLegend",
            value = p,
            pos = .GlobalEnv
        )
    }
    
    # Use saveLeg if you want to save a copy of the legend
    if (saveLeg) {
        ggsave(filename = "bivLegend.pdf", plot = p, device = "pdf",
               path = "./", width = 4, height = 4, units = "in",
               dpi = 300)
    }
    seqs <- seq(0, 100, (100 / nbreaks))
    seqs[1] <- 1
    col.matrix <- col.matrix[c(seqs), c(seqs)]
    attr(col.matrix, "breakstyle") <- breakstyle
    attr(col.matrix, "nbreaks") <- nbreaks
    return(col.matrix)
}

col.matrixQ <- colmat(nbreaks = 5, breakstyle = "quantile",
                      xlab = "Slope", ylab = "Human footprint", 
                      bottomright = "sienna4", upperright = "black",
                      bottomleft = "azure2", upperleft = "deepskyblue",
                      saveLeg = FALSE, plotLeg = TRUE)


# Function to assign colour-codes to raster data
# As before, by default assign tercile breaks
bivariate.map3 <- function(rasterx, rastery, colourmatrix = col.matrix,
                           export.colour.matrix = TRUE,
                           outname = paste0("colMatrix_rasValues", names(rasterx))) {
    # TO DO - replace raster with terra #
    require(raster)
    require(classInt)
    # export.colour.matrix will export a data.frame of rastervalues and RGB codes 
    # to the global environment outname defines the name of the data.frame
    quanx <- getValues(rasterx)
    tempx <- data.frame(quanx, quantile = rep(NA, length(quanx)))
    brks <- with(tempx, classIntervals(quanx,
                                       n = attr(colourmatrix, "nbreaks"),
                                       style = attr(colourmatrix, "breakstyle"))$brks)
    ## Add (very) small amount of noise to all but the first break
    ## https://stackoverflow.com/a/19846365/1710632
    brks[-1] <- brks[-1] + seq_along(brks[-1]) * .Machine$double.eps
    r1 <- within(tempx, quantile <- cut(quanx,
                                        breaks = brks,
                                        labels = 2:length(brks),
                                        include.lowest = TRUE))
    quantr <- data.frame(r1[, 2])
    quany <- getValues(rastery)
    tempy <- data.frame(quany, quantile = rep(NA, length(quany)))
    brksy <- with(tempy, classIntervals(quany,
                                        n = attr(colourmatrix, "nbreaks"),
                                        style = attr(colourmatrix, "breakstyle"))$brks)
    brksy[-1] <- brksy[-1] + seq_along(brksy[-1]) * .Machine$double.eps
    r2 <- within(tempy, quantile <- cut(quany,
                                        breaks = brksy,
                                        labels = 2:length(brksy),
                                        include.lowest = TRUE
    ))
    quantr2 <- data.frame(r2[, 2])
    as.numeric.factor <- function(x) {
        as.numeric(levels(x))[x]
    }
    col.matrix2 <- colourmatrix
    cn <- unique(colourmatrix)
    for (i in 1:length(col.matrix2)) {
        ifelse(is.na(col.matrix2[i]),
               col.matrix2[i] <- 1, col.matrix2[i] <- which(
                   col.matrix2[i] == cn
               )[1]
        )
    }
    # Export the colour.matrix to data.frame() in the global env
    # Can then save with write.table() and use in ArcMap/QGIS
    # Need to save the output raster as integer data-type
    if (export.colour.matrix) {
        # create a dataframe of colours corresponding to raster values
        exportCols <- as.data.frame(cbind(
            as.vector(col.matrix2), as.vector(colourmatrix),
            t(col2rgb(as.vector(colourmatrix)))
        ))
        # rename columns of data.frame()
        colnames(exportCols)[1:2] <- c("rasValue", "HEX")
        # Export to the global environment
        assign(
            x = outname,
            value = exportCols,
            pos = .GlobalEnv
        )
    }
    cols <- numeric(length(quantr[, 1]))
    for (i in 1:length(quantr[, 1])) {
        a <- as.numeric.factor(quantr[i, 1])
        b <- as.numeric.factor(quantr2[i, 1])
        cols[i] <- as.numeric(col.matrix2[b, a])
    }
    r <- rasterx
    r[1:length(r)] <- cols
    return(r)
}



# create the bivariate raster
bivmapQ <- bivariate.map3(rasterx = slp, rastery = hft,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixQ)

setwd("./output/")
writeRaster(bivmapQ, "./outputs/Bivar.tif")


library(MASS) 
library(reshape2) 
library(reshape) 
library(data.table)
library(ggpro)

bivMapDFQ <- setDT(as.data.frame(bivmapQ, xy = TRUE))
colnames(bivMapDFQ)[3] <- "BivValue"

mex
rich2
clipExt <- extent(-118.5, -86.5, 14, 33)
bivMapDFQ <- melt(bivMapDFQ, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

# Make the map using ggplot
map_q <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
    geom_raster(aes(fill = bivVal)) +
    scale_y_continuous(breaks = seq(14, 34 , by = 8), 
                       labels = paste0(seq(14, 34, 8), "°")) +
    scale_x_continuous(breaks = seq(-119, -86.5 ,8), 
                       labels = paste0(seq(-119, -86.5 ,8), "°")) +
    scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
    theme_classic() +
    theme(text = element_text(size = 10, colour = "black")) +
    borders(colour = "black", size = 0.5) +
    coord_quickmap(expand = FALSE, xlim = clipExt[1:2], ylim = clipExt[3:4]) +
    theme(legend.position = "none",
          plot.background = element_blank(),
          strip.text = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black")) +
    labs(x = " ", y = " ")


ggsave(plot = map_q,
       filename = "./outputs/BivariatePlot3.tiff",
       device = "tiff", width = 20, height = 17, units = "cm",
       dpi = 400)

ggsave(plot=col.matrixQ, filename = "colmatrix.tiff", device = "tiff", width = 5, height = 5, units = "cm", dpi=400)


# Otra forma --------------------------------------------------------------

# Mapa bivariado ----------------------------------------------------------
#ejemplo curso fabricii
library(terra)
library(sf)

slp<-rast("D:/Mauro/borrar/wc2.1_30s_slope_MX.tif")
hfp<-rast("D:/Mauro/borrar/bivar/hfp2013_merisINT_MX.tif")
riq<-rast("D:/Mauro/borrar/Riqueza_MX_HFP.tif")

hfp.res<- hfp |> 
  resample(slp, method= "near")

slp.vec <- slp |> 
  as.polygons(dissolve = F) |> 
  st_as_sf() 

hfp.vect <- hfp.res |> 
  as.points() |> 
  st_as_sf()

riqxslp<-st_join(slp.vec, hfp.vect)
names(riqxslp)<- c("Slope", "HFP", "geometry")

st_write(riqxslp, "D:/Mauro/borrar/bivar/hfpxslp.shp", append=FALSE)

d=expand.grid(x=1:3,y=1:3)
names(d)<-names(riqxslp)[c(1,2)]


leg= ggplot(d, aes(x=Slope,y=HFP))+
  geom_tile(aes(alpha=Slope+HFP,fill=atan(HFP/Slope)))+
  scale_fill_viridis_c()+
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  coord_equal()

ggsave(filename = "D:/Mauro/borrar/bivar/biv_slp_hfp_leg.png",
       width = 15,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 300)


x=classInt::classIntervals(riqxslp$Slope,4,style = "quantile")
y=classInt::classIntervals(riqxslp$HFP,4,style = "quantile")

riqxslp$x = classInt::findCols(x)
riqxslp$y = classInt::findCols(y)
riqxslp$alpha = as.character(riqxslp$x + riqxslp$y)
riqxslp$color = as.character(atan(riqxslp$y/riqxslp$x))

map = ggplot()+
  geom_sf(data = riqxslp,aes(fill=color,alpha=alpha),shape=15, size=11,show.legend = F, color = NA)+
  scale_fill_viridis_d()+
  theme_void()

ggsave(filename = "D:/Mauro/borrar/bivar/biv_slp_hfp.tiff",
       width = 20,
       height = 17, #alto
       units ="cm",
       dpi = 400)


library(patchwork)
map + leg


# Sankey diagrams ---------------------------------------------------------

#Land use-land cover change plot

# upload necessary packages
library(sf)        
library(tmap)      
library(pals)      
library(dplyr)     
library(tidyr)     
library(ggplot2)   
library(networkD3) # Sankey diagram

library(terra)
CC<-rast("./INEGI/Cambios/TC/CC_TC_ras.tif")
CT<-rast("./INEGI/Cambios/TC/CT_TC_ras.tif")
TC<-rast("./INEGI/Cambios/TC/TC_TC_ras.tif")
TT<-rast("./INEGI/Cambios/TC/TT_TC_ras.tif")

TT |> 
  values() |> 
  na.omit() |> 
  length()

links <- data.frame(
  source=c("C2","C2", "T2", "T2"), 
  target=c("C7","T7", "C7", "T7"), 
  value=c(1105003, 187494, 127372, 1054470)
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1


sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, fontSize = 25, 
              nodeWidth = 12,
              iterations=0,
              margin = list("left" = 100))


