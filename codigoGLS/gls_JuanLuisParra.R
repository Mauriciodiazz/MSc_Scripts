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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Datos
data <- read.csv("data.csv")
names(data)
data

# Normalized variables
# z-transform on our data (mean-centering and dividing by one standard deviation)
datos_z <- decostand(data[,c(4:7)], "standardize") # indicamos las variables que vamos a transformar (var explicativas) esto hace lo mismo que la funcion scale
datos_z <- cbind(data[,c(1:3,8)], datos_z) # unimos el resto de variales
head(datos_z)

# test of normality
# Shapiro-Wilk's test
shapiro.test(datos_z$y) # verificamos normalidad de la variable respuesta - No hay normalidad en la variable respuesta
format(shapiro.test(datos_z$y)$p.value, scientific = F)
hist(datos_z$y)

##################################################
#Spatial correlation test
#Moran's Test
dists <- as.matrix(dist(cbind(datos_z$coordx, datos_z$coordy))) # es importante que no existan coordenadas duplicadas! o sitios con las mismas coordenadas. Esto es una matriz de distancias cartesianas entre cada par de coordenadas que repreentan un punto en el espacio
dists.inv <- 1/dists #inverso de esa distancia
diag(dists.inv) <- 0 #la diagonal de esa matriz tiene un valor de cero porque la distancia al mismo punto es cero y este objeto tiene valores de infinito
dists.inv[1:5, 1:5]
##
moranI <- Moran.I(x=datos_z$y, weight=dists.inv) # p.value significativo indica fuerte autocorrelacion espacial
#La matriz de pesos es la distancia euclidiana que hay entre cada ímtp
moranI

#Spatial correlation correction with GLS
#function using the correlation argument. We fit our model using different correlation structures, and we then use AIC to choose the best model. 
#Note that, as with the Variogram function, we need two columns in our dataframe containing the coordinates of our sites. The nugget argument allows us to choose wether we want a nugget effect (intercept) or not.

model1 <- gls(log(y)~x1+x2+x3+x4, # transformamos "y" con log para garantizar normalidad 
               correlation = corExp(form = ~coordy+coordx, nugget = TRUE), 
               data = datos_z,
              control = glsControl(opt = "optim"))


model2 <- gls(log(y)~x1+x2+x3+x4, 
              correlation = corGaus(form = ~coordy+coordx, nugget = TRUE), 
              data = datos_z, 
              control = glsControl(opt = "optim"))

model3 <- gls(log(y)~x1+x2+x3+x4, 
              correlation = corSpher(form = ~coordy+coordx, nugget = TRUE), 
              data = datos_z,
              control = glsControl(opt = "optim"))

model4 <- gls(log(y)~x1+x2+x3+x4, 
              correlation = corLin(form = ~coordy+coordx, nugget = TRUE), 
              data = datos_z,
              control = glsControl(opt = "optim"))

model5 <- gls(log(y)~x1+x2+x3+x4, 
              correlation = corRatio(form = ~coordy+coordx, nugget = TRUE), 
              data = datos_z,
              control = glsControl(opt = "optim"))

#Check what model has best AIC
AICc(model1, model2, model3, model4, model5)
#summary best models
summary(model4)


# ---- FIN