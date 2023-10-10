
#We load the packages
library(stats)
library(glmmTMB)
library(ggplot2)
library(lme4)
library(scales)
library(nlme)
library(caret)
library(gridExtra)
library(MuMIn)
library(car)
library(multcomp)
library(gmodels)
library (spdep)
library (ncf)
library (ade4)
library (gstat)
library (geoR)
library (sp)
library (Hmisc)
library (nlme)
library(DHARMa)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

####GROUND HUNTERS#######
#Load data #Ground hunters
Datos_GH<- read.csv("Ground_hunters.csv")
Datos_GH
attach(Datos_GH)
names(Datos_GH)


#Scale landscape variables
X.forest_300mB <- scale(X.forest_300m, center = F)


# Ground hunters Abundance
ModelAbu_GH <- dredge (glmer.nb (Abundance~ Diversification + X.forest_300mB +  Stage + Diversification*X.forest_300mB + Stage*Diversification + Stage*X.forest_300mB + (1|Plots), data=Datos_GH, na.action = "na.fail"))

ModelAbu_GH
#The most probable models are selected based on the AIC. See Burnham and Anderson, 2002.


#Probable models
ModProb_Abu1_GH <- glmer.nb(Abundance~ Diversification + Stage + X.forest_300mB + (1|Plots), data=Datos_GH, na.action = "na.fail")
summary(ModProb_Abu1_GH)
AIC(ModProb_Abu1_GH)

ddNew3A <-dotplot(ranef(ModProb_Abu1_GH ,condVar=TRUE))
do.call(grid.arrange,c(ddNew3A,list(nrow=1)))


ModProb_Abu2_GH <- glmer.nb(Abundance~ Stage + X.forest_300mB + (1|Plots), data=Datos_GH, na.action = "na.fail")
summary(ModProb_Abu2_GH)
AIC(ModProb_Abu2_GH)

#Null model
ModNulo_Abu_GH <- glm.nb(Abundance~ 1, data=Datos_GH)
ModNulo_Abu_GH
AIC(ModNulo_Abu_GH)

##
AIC(ModNulo_Abu_GH)#479.3605
AIC(ModProb_Abu1_GH)#441.5231 #Selected model
AIC(ModProb_Abu2_GH)#444.8216

#¿Mantenemos el factor aleatorio?

ModProb_Abu1B_GH <- glm.nb(Abundance~ Diversification + Stage + X.forest_300mB, data=Datos_GH)
summary(ModProb_Abu1B_GH)
AIC(ModProb_Abu1B_GH)


AIC(ModProb_Abu1_GH)
AIC(ModProb_Abu1B_GH)

#Se mantiene el factor aleatorio

#Selected model
ModProb_Abu1B_GH <- glm.nb(Abundance~ Diversification + Stage + X.forest_300mB, data=Datos_GH)
summary(ModProb_Abu1B_GH)


#The model residuals are analyzed
#let's see the residuals
par(mfrow=c(2,2))
plot(ModProb_Abu1B_GH)

shapiro.test(residuals(ModProb_Abu1B_GH))

#Overdispersion assumption are evaluated
testDispersion(simulateResiduals(ModProb_Abu1B_GH))
#It`s OK

#The spatial autocorrelation of the residuals is analyzed through variography

XY<- cbind(Datos_GH$X_coord,Datos_GH$Y_coord)
XY
#Coordinates are converted to a matrix
XYcoords<- as.matrix(XY)
XYcoords
#Coordinates are centered
XYcentre <- scale(XYcoords, scale = FALSE)
XYcentre
X <- XYcentre[,1]
Y <- XYcentre[,2]
XY<-cbind(X,Y)
XY
#
XYcoords<-as.matrix(XY)
XYcoords
XYcoords<-data.frame(XYcoords)
#Repeated positions are searched
dup.coords(XYcoords)
#The value 0.001 is used to remove duplicate coordinates
XYcoords_SD <-jitterDupCoords(XYcoords[1:2], max= 0.1)
dup.coords(XYcoords_SD)#There are no more repeated coordinates

#Variograms
#The matrix is made with the coordinates and the model residuals
Variography_Mo_Ab_GH <- data.frame(XYcoords_SD,(residuals(ModProb_Abu1B_GH)))
Variography_Mo_Ab_GH

#Transforms into geodata format
Variography_Mo_Ab_GH_Geo<-as.geodata(Variography_Mo_Ab_GH)
Variography_Mo_Ab_GH_Geo

#The distance for the variogram is calculated, 60% of the distance
sort(dist(Variography_Mo_Ab_GH[1:2]))
max(dist(Variography_Mo_Ab_GH[1:2]))
porcen_60 <- (max(dist(Variography_Mo_Ab_GH[1:2]))/100)*60
porcen_60

#Variogram
#Plot semivariogram
Vgm_Mo_Ab_GH<-variog(Variography_Mo_Ab_GH_Geo,max.dist=500)
vgm_cloud_GH<-variog(Variography_Mo_Ab_GH_Geo, option="cloud", max.dist=500)
par(mfrow=c(1,2))
plot(vgm_cloud_GH, main="Abundance cloud",xlab="Distance (m)",ylab="Semi-variance")
plot(Vgm_Mo_Ab_GH,main="Abundance variogram with CI",xlab="Distance (m)",ylab="Semi-variance")

#Theoretical models fit
###Spherical model###
sphp_Ab_GH<-variofit(Vgm_Mo_Ab_GH,cov.model="sph", weights="npairs")
summary(sphp_Ab_GH)
sphp_Ab_GH$cov.pars[1] + sphp_Ab_GH$nugget #Sill
sphp_Ab_GH$nugget#nugget
wtss_Ab_GH<-sum(Vgm_Mo_Ab_GH$n*(Vgm_Mo_Ab_GH$v - weighted.mean(Vgm_Mo_Ab_GH$v,Vgm_Mo_Ab_GH$n))^2)#weighted least squares are calculated to calculate the fit of the theoretical model
Sphwr2_Ab_GH<- 1 - sphp_Ab_GH$value / wtss_Ab_GH #model fit (R^2)
Sphwr2_Ab_GH

###Exponetial model###
expp_Ab_GH <- variofit (Vgm_Mo_Ab_GH, cov.model = "exp", weights = "npairs")
summary(expp_Ab_GH)
expp_Ab_GH$cov.pars [1] + expp_Ab_GH$nugget#Sill
expp_Ab_GH$nugget#nugget
exppr2_Ab_GH<- 1 - expp_Ab_GH$value / wtss_Ab_GH #model fit (R^2)
exppr2_Ab_GH

###Gaussian model###
gaup_Ab_GH <- variofit (Vgm_Mo_Ab_GH, cov.model = "gau", weights = "npairs")
summary(gaup_Ab_GH)
gaup_Ab_GH$cov.pars [1] + gaup_Ab_GH$nugget#Sill
gaup_Ab_GH$nugget#nugget
gaupr2_Ab_GH<- 1 - gaup_Ab_GH$value / wtss_Ab_GH #model fit (R^2)
gaupr2_Ab_GH

#There is no significant patch of aggregation in the model residuals
#Finish

