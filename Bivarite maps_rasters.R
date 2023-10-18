# Mapa bivariado ----------------------------------------------------------
# voy a hacer el mapa a 10km2 no a 1km2
#ejemplo curso fabricio
library(terra)
library(sf)
library(letsR)
library(tidyverse)
library(patchwork)
library(biscale)


# Cargar capa de pendiente ------------------------------------------------
slp<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")


# Raster de riqueza -------------------------------------------------------
#cargar matriz PAM
spd.list.join<- read.table("./outputs/tablas/spp_PAM.txt",sep="\t", dec=".", header=T)
head(spd.list.join)

# Crear la PAM a partir de los puntos de presencia a una resolucion de 10km
PAM<-lets.presab.points(xy= as.matrix(spd.list.join[,1:2]),
                        species= spd.list.join$spp,
                        xmn = -120,
                        xmx = -85,
                        ymn = 14,
                        ymx = 32, #10km 0.08333 / 30km 0.24999 /50km 0.41665
                        resol = 0.08333) 

plot(PAM)
riq<-PAM$Richness_Raster
plot(riq)


# Remuestrear slope y volverla poligono ---------------------------------
slp.vec <- slp |> 
  resample(riq, method="bilinear") |> 
  as.polygons(dissolve = F) |> 
  st_as_sf() 

riq.vect <- riq |> 
  as.points() |> 
  st_as_sf()

riqxslp<-st_join(slp.vec, riq.vect)
names(riqxslp)<- c("slope", "z", "geometry")


data <- bi_class(riqxslp, x = z, y = slope, style = "quantile", dim = 3)
#data <- bi_class(stl_race_income, x = pctWhite, y = medInc, style = "quantile", dim = 3)

map <- 
  ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), colour = NA, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 3) +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


legend <- 	
  bi_legend(
    pal = "GrPink",
    dim = 3,
    xlab = "Terrain slope",
    ylab = "Species richness",
    size = 8, arrows=F)

map + inset_element(legend, 0, 0, 0.3, 0.3, align_to="full")

ggsave(filename = "./outputs/biv_slp_riq2.png",
         width = 15,
         height = 10, #alto
         scale=2,
         units ="cm",
         dpi = 300)


#mapa solo riqueza
plot(riq)
install.packages("tidyterra")
library(tidyterra)

#mapa mexico, de referencia
mexico <- st_read("./Mexico_Estados/México.shp")

mask(riq, mexico) |> 
  plot()

rich.map<-ggplot()+
  geom_spatraster(data=riq |>
                    subst(0,NA) |> 
                    mask(mexico)) +
  scale_fill_viridis(option="turbo",na.value = 'transparent')+
  geom_sf(data=mexico, fill=NA)+
  labs(fill="Species richness")+
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(), legend.position = c(0.2,0.2))

((map + inset_element(legend, 0, 0, 0.3, 0.3, align_to="full")) / rich.map)

#qué rango de riqueza tiene cada categoría?
head(data)
data |> as_tibble() |> 
  select(z,slope,bi_class) |> 
#  filter(bi_class=="1-1") |> 
  group_by(bi_class) |> 
  summarise(min.z=min(z),
            max.z=max(z),
            min.slp=min(slope),
            max.slp=max(slope))


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


