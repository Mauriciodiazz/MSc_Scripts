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
PAM.list.join<- read.table("./outputs/tablas/spp_PAM.txt",sep="\t", dec=".", header=T)
head(PAM.list.join)

# Crear la PAM a partir de los puntos de presencia a una resolucion de 10km
PAM<-lets.presab.points(xy= as.matrix(PAM.list.join[,1:2]),
                        species= PAM.list.join$spp,
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

plot(riqxslp)

data <- bi_class(riqxslp, x = z, y = slope, style = "quantile", dim = 3)
#data <- bi_class(stl_race_income, x = pctWhite, y = medInc, style = "quantile", dim = 3)
plot(data)
data$class2<-data$bi_class
data$class2[which(data$bi_class=="1-1")]<-1
data$class2[which(data$bi_class=="1-2")]<-2
data$class2[which(data$bi_class=="1-3")]<-3
data$class2[which(data$bi_class=="2-1")]<-4
data$class2[which(data$bi_class=="2-2")]<-5
data$class2[which(data$bi_class=="2-3")]<-6
data$class2[which(data$bi_class=="3-1")]<-7
data$class2[which(data$bi_class=="3-2")]<-8
data$class2[which(data$bi_class=="3-3")]<-9

data2<-vect(data)
writeVector(data2,"./outputs/temps/by_class.shp")
riq2<-riq
writeRaster(riq2,"./outputs/temps/riq2.tif")
a<-rasterize(data2, riq2, values=data2$class2)
a


a<-vect(data)
a2<-terra::extract(a, riq.vect)


### DATA: 
### Delta = a presence-absence matrix (PAM); ***species are rows and sites are columns*** (without coordinates or any other information)
### *Note that most PAMs nowadays come as transpose; with species as columns and sites as rows. 
### *If so, you'd need to transpose your PAM before using this function. 


PAM<-matrix(data=c(0,0,1,0,1,0,1,1,1), ncol = 3, nrow=3)
row.names(PAM)<-c("sp1","sp2","sp3")
colnames(PAM)<-c("s1","s2","s3")
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

#### by SITES ###

sites.RDplot<- function(pamRD) {
  
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
  plot(SiteRange/Quadrats, SpeciesRichness/Species, xlim = c(0,1), ylim = c(0,1), xlab = "Mean proportional per-site range", ylab = "Proportional species richness", pch=21,bg="black",cex=1.5)
  
  
  
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



RangeDiversity(PAM) |> 
  sites.RDplot()

sum(PAM)
1/5/9



library(terra)
library(tidyverse)
area<-rast("./WorldClim_30s/wc2.1_30s_elev/Area_slp_cor_AG/AC_buf.tif")
slp<-rast("WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")

#1. reproyectar capa de area
areawgs<-project(area, slp, method="bilinear")

#2. extraer valores de pendiente y area
slp.area<-extract(areawgs, as.points(slp), bind=T) |> 
  as.data.frame()

slp.area |> 
  na.omit() |> 
 # sample_n(size=100) |> 
  mutate(AC_buf=AC_buf/1000000, #a km cuadrado
         slope_clas=cut(slope_mx_g , breaks = 0:31, labels = 0:30)) |> 
  group_by(slope_clas) |> 
  summarise(sum=sum(AC_buf)) |>
  na.omit() |> 
  ggplot(aes(x=slope_clas, y=sum))+
  geom_bar(stat='identity')+
  labs(x="Terrain slope", y="Total area (km\u00b2)") +
  theme_classic()


1000*1000 / 1000000

1000000






library(terra)
point<-vect("./borrar/ablancae_sites.shp")
plot(point)
distance(point) |> 
  mean()
which(dist>29000)

library(tidyverse)
data<-read.table("C:/Users/Mauricio Diaz/Documents/Pregrado/Tesis_preg/Muestreos/especies_sitios_curada.txt", sep="\t", dec=".", header=T)
head(data)

data |> 
  filter(names=="20210419_062000.WAV" & sitio== "26")


tocopy<-data |> 
  mutate(coun=1) |> 
  group_by(names, sitio) |> 
  summarise(num.spp=sum(coun)) |> 
  arrange(desc(num.spp)) |> 
  mutate(names=paste0(sitio,"_",names))

tocopy<-tocopy[1:80,]
write.table(tocopy,"C:/Users/Mauricio Diaz/Documents/Trabajos/Proyectos/birdnet/audios_birdnet/spp_audios.txt",sep="\t", dec=".", row.names = F)

a<-data[data$names %in% tocopy$names,] 
write.table(a,"C:/Users/Mauricio Diaz/Documents/Trabajos/Proyectos/birdnet/audios/spp_sites.txt", sep="\t", dec=".", row.names = F)

a |> 
  filter(names=="20210422_062000.WAV" & sitio=="68")

data[which(data$names=="20210515_054800.WAV"),]

my.file.copy <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.copy(from = from,  to = to)
}


for (x in 1:nrow(tocopy)) {
  my.file.copy(from = paste0("G:/Grabaciones_Ablancae/","Sitio ",tocopy[x,2],"/","Tesis ",tocopy[x,2],"/",tocopy[x,1]),
               to = paste0("C:/Users/Mauricio Diaz/Documents/Trabajos/Proyectos/birdnet/audios/",tocopy[x,2],"_",tocopy[x,1]))
}

tocopy |> 
  filter(sitio==1)


mx<-vect("./Mexico_Estados/Mexico_continent.shp")
slp1k<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")

slop250<-rast("C:/Users/Mauricio Diaz/Downloads/borrar/slope_90M_n00w090/dtm_slope.tif")
roug<-rast("C:/Users/Mauricio Diaz/Downloads/borrar/slope_90M_n00w090/dtm_roughness.tif")

slop2<- slop250 |> 
  crop(mx) |> 
  mask(mx) |> 
  resample(slp1k)

roug2<- roug |> 
  crop(mx) |> 
  mask(mx) |> 
  resample(slp1k)

plot(a)

a<-c(slop2, roug2, slp1k) |> 
  as.data.frame() |> 
  na.omit()
head(a)

cor(a, method="spearman", use = "na.or.complete")
cor(a, method="spearman")

a |> 
  ggplot(aes(x=slope_mx_g , y=dtm_roughness))+
  geom_point()+
  labs(x="Terrain slope", y="Terrain roughness")+
  theme_classic()

z<-rast("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/Riqueza/Riqueza_MX.tif")
inegi<-rast("./INEGI/Cambios/INEGI_VII_TC.tif") |> 
  resample(z, method="near") |> 
  subst(2,100)
  
inegi<-subst(inegi, 2, 100)

a<-z+inegi
plot(a)
writeRaster(a, "./borrar/Zxcat.tif")


