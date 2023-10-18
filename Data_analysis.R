#Analisis riqueza vs pendiente

Z.slp.cats<- read.table("./outputs/Z.slp.cats.txt", sep = "\t", dec=".", header=T)

library(lme4)
library(MASS)
library(moments)
library(tidyverse)
library(easystats)

#THETA: PARAMETRO DE LA VARIACIÓN DE LAS OBSERVACIONES; ESTIMA QUE TAN GRANDE ES LA VARIANZA RESPECTO AL PROMEDIO
Z.slp.cats$cat<-factor(Z.slp.cats$cat)
Z.slp.cats$ecoreg<-factor(Z.slp.cats$ecoreg)

options(contrasts  =c("contr.treatment", "contr.poly"))
mblrdata.dist <- rms::datadist(Z.slp.cats)
options(datadist = "Z.slp.cats")

mod.null= glm(z ~ slope*cat,
              data = Z.slp.cats,
              family = neg.bin(theta=log(var(Z.slp.cats$z)/mean(Z.slp.cats$z)))) 

mod = glmer(z ~ slope*cat + (1|ecoreg), 
            data = Z.slp.cats,
            family = neg.bin(theta=log(var(Z.slp.cats$z)/mean(Z.slp.cats$z))))

mod2 = glmer(z ~ slope*cat + (1 + cat | ecoreg), 
            data = Z.slp.cats,
            family = neg.bin(theta=log(var(Z.slp.cats$z)/mean(Z.slp.cats$z))))

mod.qp = glmmPQL(z ~ slope*cat, random= ~ 1|ecoreg, 
            data = Z.slp.cats,
            family = quasipoisson(link='log'))

mod2.qp = glmmPQL(z ~ slope*cat, random=  ~ 1 + cat | ecoreg, 
            data = Z.slp.cats,
            family = quasipoisson(link='log'))

# Reporte del modelo ------------------------------------------------------
report(mod)

#¿cuáles implicaciones tiene usar pendientes aleatorias pra lo que quiero hacer?
#¿cual es la diferencia entre quasipoisson y binomial negativa? 
#Por qué QP no tiene AIC?
#Por qué QP no saca pseudo R2?


# sobredisperson entre los modelos ----------------------------------------

# quasi-poison ------------------------------------------------------------
# extract pearson residuals
PearsonResiduals <- resid(mod2.qp, type = "pearson")
# extract number of cases in model
Cases <- nrow(Z.slp.cats)
# extract number of predictors (plus intercept)
NumberOfPredictors <- length(fixef(mod2.qp)) +1
# calculate overdispersion
Overdispersion.qp2 <- sum(PearsonResiduals^2) / (Cases-NumberOfPredictors)
# inspect overdispersion

Overdispersion.qp
Overdispersion.qp2
Overdispersion.bn
Overdispersion.bn2

performance::check_overdispersion(mod)
r2.mod<-r2_nakagawa(mod)
# Conditional R2: 0.587
# Marginal R2: 0.017
summary(mod)
r2.mod2<-r2_nakagawa(mod2) 
# Conditional R2: 0.576
# Marginal R2: 0.018

r2_nakagawa(mod.qp)
r2_nakagawa(mod2.qp)
performance::r2(mod.qp)

# negative binomial -------------------------------------------------------
# extract pearson residuals
PearsonResiduals <- resid(mod, type = "pearson")
# extract number of betas + predictors + sigma
NumberOfPredictors <- 2+1+1
# extract number of cases in model
Cases <- nrow(Z.slp.cats)
# calculate overdispersion parameter
Overdispersion.bn <- sum(PearsonResiduals^2) / (Cases / NumberOfPredictors)# show overdispersion parameter
Overdispersion.bn
Overdispersion.bn2


AIC(logLik(mod.null))
AIC(logLik(mod))
AIC(logLik(mod.qp))
AIC(logLik(mod2)) #mejor modelo!

# test random effects
null.id = -2 * logLik(mod.null) + 2 * logLik(mod)
pchisq(as.numeric(null.id), df=1, lower.tail=F) 

anova(mod, mod2, test="Chi")

# plot model
plot_model(mod2, type = "pred", terms = c("slope", "cat"))

summary(mod)

mod.null = glmer(z ~ 1 + (1|ecoreg), 
            data = Z.slp.cats,
            family = neg.bin(theta=log(var(Z.slp.cats$z)/mean(Z.slp.cats$z))))

mb.lmer = lmer(z ~ slope*cat + (1 + cat | ecoreg), REML = T, data = Z.slp.cats)
plot(mod2, ecoreg ~ resid(.), abline = 0 )

# test vifs
car::vif(mb.lmer)
car::vif(mod)

1-var()
performance::check_overdispersion(mod2) #Overdispersion detected. 
performance::check_overdispersion(mb.lmer) #Overdispersion detected. 

resid<- resid(mod)
res.mod<-mod |> 
  resid() |> 
  as_tibble() |> 
  #select(value) |> 
  pull()

res.modnull<-mod.null |> 
  resid() |> 
  as_tibble() |> 
 # select(value) |> 
  pull()

#pseudo R cuadrado
1-(var(res.mod)/var(res.modnull))

r2_nakagawa(mod)
# Conditional R2: 0.587
# Marginal R2: 0.017

install.packages("MuMIn")
MuMIn::r.squaredGLMM(mod, null=mod.null)

library(easystats)
report(mod)

save.image()
#vectores para las pendientes de las relaciones dadas por el glmer. Lo que pasa es que las pendientes asociadas al glmer para cada categoría estan en relación de la primera (CC), es decir que hay que re calcular la pendiente, entonces se hace el siguiente paso


performance::check_overdispersion(mb.lmer)

CC.sl<- c(0,1,0,0,0,0,0,0)
CT.sl<- c(0,1,0,0,0,1,0,0)
TC.sl<- c(0,1,0,0,0,0,1,0)
TT.sl<- c(0,1,0,0,0,0,0,1)

#luego se realiza esta prueba, en donde se comparan las pendientes por pares de categorias. Es una t de student
#t= resta de pendientes/ sde -> valor de distribucion de probabilidad
library(multcomp)
#HO:No hay diferencias entre las pendientes
summary(glht(mod,rbind(CC.sl-CT.sl, 
                       CC.sl-TC.sl, 
                       CC.sl-TT.sl,
                       CT.sl-TC.sl, 
                       CT.sl-TT.sl,
                       TC.sl-TT.sl,
                       CC.sl,
                       CT.sl,
                       TC.sl,
                       TT.sl))) |> 
  broom::tidy() |> 
  mutate(contrast=c("CC vs CT", 
                    "CC vs TC", 
                    "CC vs TT",
                    "CT vs TC", 
                    "CT vs TT",
                    "TC vs TT",
                    "CC pendiente",
                    "CT pendiente",
                    "TC pendiente",
                    "TT pendiente"),
         signif = gtools::stars.pval(adj.p.value)) |> 
  kableExtra::kbl(digits = 3) |> 
  kableExtra::kable_styling()


# ACUALIZACIÓN! -----------------------------------------------------------

#Juan me hizo caer en cuenta de que es probable que exista 1. autocorrelación espacial y que además, la pendiente esté correlacionada con las ecoregiones 2. Hay una forma de solucionarlo y es aplicando GLS o usando eigenvectors spatial filters....

#1. Voy a ver como se comporta la pendiente dentro de cada ecoregión a ver si hay o no correlación. Sin embargo, eso no elimina la autocorrelación espacial, por lo que luego evaluaré los GLS y si no, pues los eigenvectors

head(Z.slp.cats)
Z.slp.cats |>
#  sample_n(size=100) |> 
  ggplot(aes(x=as.character(ecoreg), y=slope))+
  geom_boxplot()
  

levels(factor(Z.slp.cats$ecoreg)) |> 
  length()
#Si hay una relación entre la pendiente y las ecoregiones

# problema: Para hacer los GLS y evaluar autocorrelación espacial necesito crear una matriz de distancias entre los pixeles, además, parar ealizar el índice de moran es necesario hacerlas. Por lo tanto, para facilidad en los análisis voy a reescalar el mapa de riqueza a 5km.


# Mapa de riqueza con LetsR -----------------------------------------------
#Direcciones de los rasters de distribución potencial
ras.dir<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = T, pattern =".tif$")
#Lista de nombres de los archivos
ras.list<-list.files("F:/Maestria_DD/spp.records_DD/specialists_DD/spec.shapes_DD/Modelos_binarios/MB_Mex/", full.names = F, pattern =".tif$") |> 
  substr(1,7) #Esto solo elimina los caracteres .tif


# Crear PAM y raster de riqueza con Lets R a partir de rasters ------------
library(tidyverse)
library(letsR) # lets.presab.points

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
# write.table(spd.list.join,"F:/Maestria_DD/spp.records_DD/specialists_DD/temporales/spp_PAM.txt", sep="\t", dec=".")
head(spd.list.join)

# Cargo la tabla para optimizar tiempo
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

x11()
plot(PAM)
summary(PAM)
PAM.raster<-PAM$Richness_Raster

# Cargar y remuestrear el raster de pendiente -------------------------------

slope<-rast("./WorldClim_30s/wc2.1_30s_elev/slope_mx_g_res.tif")
slope.3<-resample(slope, PAM.raster, method="bilinear")


# Cargar y remuestrear las áreas conservadas y transformadas --------------

inegivii<-rast("./INEGI/Cambios/INEGI_VII_TC.tif")
inegivii.3<-resample(inegivii, PAM.raster, method="near")


# Extraer pendiente e inegi vii para cada pix de riqueza --------
z.slp<-extract(c(slope.3,inegivii.3, PAM.raster), 
               as.points(PAM.raster), xy=T) |> 
  na.omit()

# ordenar
z.slp<-z.slp[,c(5,6,4,2,3)]
names(z.slp)<-c("x","y","z","slope","cat")
names(z.slp)
z.slp$cat<-as.character(z.slp$cat)
str(z.slp)
levels(factor(z.slp$cat))

# Evaluando estructura espacial -------------------------------------------

# semivariograma
distancias<-dist(z.slp[,c(1,2)])
summary(distancias)
head(z.slp)

##################################################
### script de Spatial Filters curso de macro

#Agregando el “espacio” al análisis de la riqueza

#Cargar los siguientes paquetes:
library(letsR)
library(sf)
library(terra)
library(ape)
library(vegan)

#Ajustar el modelo global (OLS) entre la riqueza de aves y las variables ambientales
aves.ras1.lm <- lm(z ~ slope*cat, data = z.slp)
summary(aves.ras1.lm)

#Checar la autocorrelación espacial en los residuos del modelo

#Necesitamos crear una matriz de distancias geográficas entre los sitios/celdas
aves.ras1.coords.dist <- lets.distmat(as.matrix(z.slp[, c(1, 2)]))
# Ahora sí podemos calcular el I de Moran para diferentes clases de distancia y ver el correlograma resultante
aves.ras1.lm.moran <- lets.correl(residuals(aves.ras1.lm), aves.ras1.coords.dist, 10)

# hay estructura espacial en los residuales!

#Y hay autocorrelación espacial dentro de mi variable?

# Spatial correlation test - Moran's Index --------------------------------
library(gstat)
library(lattice)
library(nlme)
library(MuMIn)
library(ape)
library(vegan)

#Moran's Test
dists <- as.matrix(dist(cbind(z.slp$x, z.slp$y))) # es importante que no existan coordenadas duplicadas! o sitios con las mismas coordenadas. Esto es una matriz de distancias cartesianas entre cada par de coordenadas que repreentan un punto en el espacio
dists.inv <- 1/dists #inverso de esa distancia
diag(dists.inv) <- 0 #la diagonal de esa matriz tiene un valor de cero porque la distancia al mismo punto es cero y este objeto tiene valores de infinito
dists.inv[1:5, 1:5]
##
moranI <- Moran.I(x=z.slp$z, weight=dists.inv) # p.value significativo indica fuerte autocorrelacion espacial
#La matriz de pesos es la distancia euclidiana que hay entre cada ímtp
moranI

#hay normalidad?
library(tidyverse)
# test of normality
# Shapiro-Wilk's test
format(shapiro.test(z.slp|>  # verificamos normalidad de la variable respuesta
                      sample_n(size=5000) |> 
                      select(z) |> 
                      pull())$p.value, scientific = F)

hist(z.slp$z) # No hay normalidad en la variable respuesta


#ajustar modelos gls
#####################################################
#modelos gls

#Load packages
library(gstat)
library(lattice)
library(nlme)
library(MuMIn)
library(ape)
library(vegan)

# Direstorio de trabajo
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Normalized variables
# z-transform on our data (mean-centering and dividing by one standard deviation) como tengo solo una variable continua, no es necesario estandarizar
# datos_z <- decostand(data[,c(4:7)], "standardize") # indicamos las variables que vamos a transformar (var explicativas) esto hace lo mismo que la funcion scale
# datos_z <- cbind(data[,c(1:3,8)], datos_z) # unimos el resto de variales
# head(datos_z)

# Ajustar un gls
#Spatial correlation correction with GLS
#function using the correlation argument. We fit our model using different correlation structures, and we then use AIC to choose the best model. 
#Note that, as with the Variogram function, we need two columns in our dataframe containing the coordinates of our sites. The nugget argument allows us to choose wether we want a nugget effect (intercept) or not.

# transformamos "z" con log para garantizar normalidad 

names(z.slp)
model1.lm <- gls(log(z+1)~slope*cat,data = z.slp)

which(is.na(z.slp))
summary(z.slp)
z.slp$cat2<-z.slp$cat
z.slp[which(z.slp$cat2==1),6]<-"C"
z.slp[which(z.slp$cat2==2),6]<-"T"

model1 <- gls(log(z+1)~slope*cat, 
              correlation = corExp(form = ~y+x, nugget = TRUE), 
              data = z.slp,
              control = glsControl(opt = "optim"))


model2 <- gls(log(z+1)~slope*cat, 
              correlation = corGaus(form = ~y+x, nugget = TRUE), 
              data = z.slp, 
              control = glsControl(opt = "optim"))

model3 <- gls(log(z+1)~slope*cat, 
              correlation = corSpher(form = ~y+x, nugget = TRUE), 
              data = z.slp,
              control = glsControl(opt = "optim"))

model4 <- gls(log(z+1)~slope*cat, 
              correlation = corLin(form = ~y+x, nugget = TRUE), 
              data = z.slp,
              control = glsControl(opt = "optim"))

model5 <- gls(log(z+1)~slope*cat, 
              correlation = corRatio(form = ~y+x, nugget = TRUE), 
              data = z.slp,
              control = glsControl(opt = "optim"))

