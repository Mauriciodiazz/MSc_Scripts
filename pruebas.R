library(terra)
library(tidyverse)

slope<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
elev<-rast("./WorldClim_30s/wc2.1_30s_elev/wc2.1_30s_elev_MX.tif")
HFP<-rast("./HFP/Dryad_2020/Human_footprint_maps/hfp2013_merisINT_MX.tif")
rich<-rast("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif")
cons<-vect("./INEGI/Cambios/TC_s7/C_disuelto_s7.shp")

plot(slope)
plot(elev)
plot(rich)
plot(cons[1])

# Correlaciones -----------------------------------------------------------


ele.point<-elev |> 
  as.points()
slp.ele.pts<-terra::extract(slope, ele.point, ID=F, bind=T)
sl.el.hfp<-terra::extract(HFP, slp.ele.pts, ID=F, bind=T)
sl.el.hfp.z<-terra::extract(rich, sl.el.hfp, ID=F, bind=T)


# Slope vs elevacion ------------------------------------------------------

# National scale
slope.cor<-
  sl.el.hfp.z |> 
  as_tibble() |> 
  select(slope_mx_g,wc2.1_30s_elev_MX) |> 
  na.omit() 
cor(slope.cor$slope_mx_g, slope.cor$wc2.1_30s_elev_MX, method="spearman")
cor.test(slope.cor$slope_mx_g, slope.cor$wc2.1_30s_elev_MX, method="spearman")

p.vals.slp.total<-numeric()
for(i in 1:500){
  slp.sampled<-slope.cor[sample(1:nrow(slope.cor), 100, replace=T),] |> 
    scale() |> 
    as_tibble()
  ct<-cor.test(slp.sampled$slope_mx_g, slp.sampled$wc2.1_30s_elev_MX, method="spearman")
  p.vals.slp.total[i]<-ct$p.value
}

poolr::bonferroni(p.vals.slp.total)

# Biome scale
library(sf)
biomas<-vect("Provincias biogeograficas/Ecoregions2017_dinerstein/Biomas_2017_Dinerstein.shp")[-7,]

data.biom<-data.frame(biome=NA, cor=NA, p.val=NA)

for (x in 1:nrow(biomas)) {
data.biom[x,1]<-biomas$BIOME_NAME[x]

bio.el<-crop(elev, biomas[x,], mask=T)|> 
  as.points()

slp.ele.pts.b<-terra::extract(slope, bio.el, ID=F, bind=T)
sl.el.hfp.b<-terra::extract(HFP, slp.ele.pts.b, ID=F, bind=T)
sl.el.hfp.z.b<-terra::extract(rich, sl.el.hfp.b, ID=F, bind=T)

b.sl.el<-
  sl.el.hfp.z.b |> 
  as_tibble() |> 
  select(slope_mx_g,wc2.1_30s_elev_MX) |> 
  na.omit() 

data.biom[x,2]<-cor(b.sl.el$slope_mx_g, b.sl.el$wc2.1_30s_elev_MX, method="spearman")

p.vals.slp.b<-numeric()
for(i in 1:500){
  slp.sampled<-b.sl.el[sample(1:nrow(b.sl.el), 100, replace=T),] |> 
    scale() |> 
    as_tibble()
  ct<-cor.test(slp.sampled$slope_mx_g, slp.sampled$wc2.1_30s_elev_MX, method="spearman")
  p.vals.slp.b[i]<-ct$p.value
}

data.biom[x,3]<-poolr::bonferroni(p.vals.slp.b)$p

}


# Crop richness in conserved area -----------------------------------------

elev.crop<-
  crop(elev, cons, mask=T)

ele.point2<-elev.crop |> 
  as.points()

slp.ele.pts2<-terra::extract(slope, ele.point2, ID=F, bind=T)
sl.el.hfp2<-terra::extract(HFP, slp.ele.pts2, ID=F, bind=T)
sl.el.hfp.z2<-terra::extract(rich, sl.el.hfp2, ID=F, bind=T)

rich.cor2<-
  sl.el.hfp.z2 |> 
  as_tibble() |> 
  select(sum,wc2.1_30s_elev_MX) |> 
  na.omit()
cor(rich.cor2$wc2.1_30s_elev_MX, rich.cor2$sum, method="spearman")
cor.test(rich.cor2$wc2.1_30s_elev_MX, rich.cor2$sum, method="spearman")

p.vals.rich<-numeric()
for(i in 1:500){
  rich.sampled<-rich.cor2[sample(1:nrow(rich.cor2), 100, replace=T),] |> 
    scale() |> 
    as_tibble()
  ct<-cor.test(rich.sampled$sum, rich.sampled$wc2.1_30s_elev_MX, method="spearman")
  p.vals.rich[i]<-ct$p.value
}

poolr::bonferroni(p.vals.rich)



a<-
  sl.el.hfp |> 
  as_tibble() |> 
  na.omit() |> 
  rename(elev="wc2.1_30s_elev_MX",
         slope="slope_mx_g",
         hfp="hfp2013_merisINT_MX") |> 
  mutate( 
    slope.cod = factor(case_when(
    slope >= 0 & slope < 2 ~ "Flat",
    slope > 2 & slope < 10 ~ "Slight",
    slope > 10 & slope < 20 ~ "Moderate",
    slope > 20  ~ "High"), levels = c("Flat", "Slight", "Moderate", "High")),
    hfp.cod= factor(case_when(
      hfp == 0 ~ "no pressure",
      hfp >= 1 & hfp <= 2 ~ "Low",
      hfp >= 3 & hfp <= 5 ~ "Moderate",
      hfp >= 6 & hfp <= 11  ~ "High",
      hfp >= 12 & hfp <= 50  ~ "Very high"), levels = c("no pressure", "Low", "Moderate", "High","Very high")))

el.sl<-a |> 
  ggplot(aes(x=slope.cod, y=elev)) +
  geom_violin(fill="gray") +
  labs(x="Terrain slope", y="Elevation (m)") +
  theme_classic()

el.hfp<-a |> 
  ggplot(aes(x=hfp, y=elev)) +
  geom_point() +
  labs(x="Human footprint", y="Elevation (m)") +
  theme_classic()

el.hfp2<-a |> 
  ggplot(aes(x=hfp.cod, y=elev)) +
  geom_violin(fill="gray") +
  labs(x="Human footprint", y="Elevation (m)") +
  theme_classic()


#library(patchwork)

el.sl / el.hfp2 / hfp.2k + plot_annotation(tag_levels = "A")

ggsave(filename = "./outputs/elevation.png",
       width = 10,
       height = 6, #alto
       scale=1.5,
       units ="cm",
       dpi = 300)
# hagamos un anova con los datos que estan por encima de 2000 msnm
hfp.2k<-
  a |> 
  filter(elev>2000) |> 
    ggplot(aes(x=hfp.cod, y=elev))+
    geom_boxplot(fill="lightgray") +
    labs(x="Human footprint", y="Elevation (m)")+
    theme_classic()
#slp.2k<-
  a |> 
  filter(elev>2000) |> 
    ggplot(aes(x=slope.cod, y=elev))+
    geom_boxplot(fill="lightgray") +
    labs(x="Terrain slope", y="Elevation (m)")+
    theme_classic()

  ggsave(filename = "./outputs/elevation3.png",
         width = 10,
         height = 10, #alto
         scale=2,
         units ="cm",
         dpi = 300)
  
hist(hfp.2k$hfp)

kruskal.test(elev~hfp.cod, data=hfp.2k)


plot(a$elev, a$slope)

# hacer el bootsrap por bioma


library(sf)
library(terra)
library(tidyverse)

bufs<-read_sf("C:/Users/Mauricio Diaz/Documents/Maestria/Tesis Msc/publicacion_Msc/Submission/response/Temp/buffers.shp")
bio1<-rast("./WorldClim_30s/wc2.1_30s_bio/wc2.1_30s_bio_Ame/wc2.1_30s_bio_1_Ame.tif")

a<-st_intersection(bio1 |> 
                     as.points() |> 
                     st_as_sf(), bufs)

plot(mx)

yuc<- 
  mx[c(14,15,17),] |> st_union()

crop(slope, vect(yuc), mask=T) |> 
  as.data.frame() |> 
  as_tibble() |> 
  summarise(median=median(slope_mx_g), sd=sd(slope_mx_g))

ar_buf<-read_csv("C:/Users/Mauricio Diaz/Documents/Maestria/Tesis Msc/publicacion_Msc/Submission/response/Temp/area_buf.csv")
ar_buf |> 
  summarise(mean=mean(area2), .by = c(Polygon, area)) |> 
  mutate(pro=mean/area)

