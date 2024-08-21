#SARs

library(spatialreg)
library(sf)
library(spdep)
library(ncf)

#La idea de usar los SAR es que se incorpore la autocorrelación espacial en el modelo a través de una lista de pesos que se le agrega al modelo en un argumento. Básicamente, esto corrige la autocorrelación espacial del modelo a través de decirle que hay autocorrelación entre los pixeles vecinos.  

#Este método utiliza 3 tipos de pesos espaciales: c("W", "C", "S") y dos distancias (min y max), por lo tanto, hay que establecer 6 modelos Wmax, Wmin, Cmax, Cmin, etc. Y elegir cuál de ellos es el mejor... Cómo? Aún no sé, pero apenas llegue de campo lo resuelvo...


# Cargar datos ------------------------------------------------------------
z.slp<- read.table("./outputs/tablas/z.slp_05.txt", header=T, dec=".", sep="\t")
head(z.slp)
nrow(z.slp)


# Spacialize points -------------------------------------------------------
sp <- st_as_sf( x = z.slp, 
                coords = c("x","y"),
                crs = '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')


# Create neighbours list of class nb --------------------------------------
k1 = knn2nb(knearneigh(sp, k = 1)) 


# Calculate distances ----------------------------------------------------
dist <- unlist(nbdists(k1, sp))
max1 <- max(dist)
min1 <- min(dist)

# create a series of neighbour matrices based on different distances (distances are in km)

d_min <- dnearneigh(sp, longlat = F, d1 = 0, d2 = min1)
d_max <- dnearneigh(sp, longlat = F, d1 = 0, d2 = max1)


# Spatial weigths creation ------------------------------------------------
scheme <- c("W", "C", "S")

# Esta función coje las distancias de los vecinos y establece los pesos segun varios estilos "W", "C", o "S" (scheme[i])
# W is row standardised (sums over all links to n)
# C is globally standardised (sums over all links to n)
# S is the variance-stabilizing coding scheme

# for min distance
#spatial_weights = sw
sw_min_w <- nb2listw(d_min, zero.policy = TRUE, style = "W") 
sw_min_s <- nb2listw(d_min, zero.policy = TRUE, style = "S") 
sw_min_c <- nb2listw(d_min, zero.policy = TRUE, style = "C") 
# for max distance
sw_max_w <- nb2listw(d_max, zero.policy = TRUE, style = "W")
sw_max_s <- nb2listw(d_max, zero.policy = TRUE, style = "C")
sw_max_c <- nb2listw(d_max, zero.policy = TRUE, style = "S")


#Por si acaso la I de Moran
moran(z.slp$z,  sw_min_w, length(k1), Szero(sw_min_w))
moran(z.slp$z,  sw_min_s, length(k1), Szero(sw_min_s))
moran(z.slp$z,  sw_min_c, length(k1), Szero(sw_min_c))
moran(z.slp$z,  sw_max_w, length(k1), Szero(sw_max_w))
moran(z.slp$z,  sw_max_s, length(k1), Szero(sw_max_w))
moran(z.slp$z,  sw_max_c, length(k1), Szero(sw_max_w))


lm.simple<- lm(z ~ slope*cat2,
               data = z.slp)

cor.lm.res <- correlog(z.slp$x, z.slp$y,
                       z = residuals(lm.simple),
                       na.rm = TRUE,
                       increment = 1,
                       resamp = 1)
plot(cor.lm.res)

# Existe autocorrelación espacial


# SAR - Simultaneous Autoregressive Model ---------------------------------

error_min_w <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_min_w,
  tol = 1e-12, 
  zero.policy = T)

error_min_s <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_min_s,
  tol = 1e-12, 
  zero.policy = T)

error_min_c <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_min_c,
  tol = 1e-12, 
  zero.policy = T)

error_max_w <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_max_w,
  tol = 1e-12, 
  zero.policy = T)

error_max_s <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_max_s,
  tol = 1e-12, 
  zero.policy = T)

error_max_c <- spatialreg::errorsarlm(
  z ~ slope*cat2,
  data = z.slp,
  listw = sw_max_c,
  tol = 1e-12, 
  zero.policy = T)

#Summary  
error_min_w_s <-summary(error_min_w, Nagelkerke=TRUE)
error_min_s_s <-summary(error_min_s, Nagelkerke=TRUE)
error_min_c_s <-summary(error_min_c, Nagelkerke=TRUE)

error_max_w_s <-summary(error_max_w, Nagelkerke=TRUE)
error_max_s_s <-summary(error_max_s, Nagelkerke=TRUE)
error_max_c_s <-summary(error_max_c, Nagelkerke=TRUE)


error_max_w$coef_lm.model #Coeficientes de un modelo lineal simple
error_max_w$coefficients #Coeficientes del SAR

IC <- c(1,0,0,0,0)
IT <- c(1,0,1,0,0)

PC <- c(0,1,0,0,0)
PT <- c(0,1,0,1,0)

library(gmodels)

estimable(error_max_w, rbind(PC, PT))



AIC(error_min_w,
    error_min_s,
    error_min_c,
    error_max_w,
    error_max_s,
    error_max_c)

schemes <- rep(c(paste(scheme[i])), times = 3)
mediandef <- rep(paste0(colnames(tablelm)[2], times = 3))
slopet <- rep(paste0("slope", 1), times = 3)
models <- c("lm_mod", "error_d1", "error_d2")

results <- 
  data.frame(
    type=rep(c("min","max"), each = 3),
    style=NA,
    NK_Rsquared = NA,
    lamnda = NA,
    LR_ratio=NA,
    LR_pval=NA,
    Ase=NA, #Asymptotic standard error
    Ase_pval=NA,
    Ase_z=NA,
    Wald=NA,
    Wald_pval=NA,
    log_ll=NA,
    sigma_sq=NA,
    AIC = NA,
    AIC_lm= NA,
    Imoran=NA)

for (i in 1:length(scheme)) {
  # Spatial weights
  sw_min <- nb2listw(d_min, zero.policy = TRUE, style = scheme[i]) 
  sw_max <- nb2listw(d_max, zero.policy = TRUE, style = scheme[i]) 
  
  # minimum and maximum SAR
  error_min <- spatialreg::errorsarlm(
    z ~ slope*cat2,
    data = z.slp,
    listw = sw_min,
    tol = 1e-12, 
    zero.policy = T)
  
  error_max <- spatialreg::errorsarlm(
    z ~ slope*cat2,
    data = z.slp,
    listw = sw_max,
    tol = 1e-12, 
    zero.policy = T)
  
  # Model's summary
  error_min_s <-summary(error_min, Nagelkerke=TRUE)
  error_max_s <-summary(error_max, Nagelkerke=TRUE)
  
  results$style[i]<-scheme[i]
  results$style[i+3]<-scheme[i]
  
  results$NK_Rsquared[i]<-error_min_s$NK
  results$NK_Rsquared[i+3]<-error_max_s$NK
  
  results$lamnda[i]<-error_min_s$lambda[[1]]
  results$lamnda[i+3]<-error_max_s$lambda[[1]]
  
  results$LR_ratio[i]<-error_min_s$LR1$statistic[1]
  results$LR_ratio[i+3]<-error_max_s$LR1$statistic[1]
  
  results$LR_pval[i]<-error_min_s$LR1$p.value
  results$LR_pval[i+3]<-error_max_s$LR1$p.value

  results$Wald[i]<-error_min_s$Wald1$statistic
  results$Wald[i+3]<-error_max_s$Wald1$statistic
  
  results$Wald_pval[i]<-error_min_s$Wald1$p.value
  results$Wald_pval[i+3]<-error_max_s$Wald1$p.value
  
  results$log_ll[i]<-error_min_s$LL
  results$log_ll[i+3]<-error_max_s$LL

  results$sigma_sq[i]<-error_min_s$s2
  results$sigma_sq[i+3]<-error_max_s$s2

  results$AIC[i]<- AIC(error_min)
  results$AIC[i+3]<- AIC(error_max)

  results$AIC_lm[i]<- error_min_s$AIC_lm.model
  results$AIC_lm[i+3]<- error_max_s$AIC_lm.model
  
  results$Imoran[i]<- moran(z.slp$z,  sw_min, length(k1), Szero(sw_min))$I
  results$Imoran[i+3]<- moran(z.slp$z,  sw_max, length(k1), Szero(sw_max))$I

}

results |> 
  arrange(AIC)

# Correlogramas 
par(mfrow=c(1,2))
cor.sar1.res <- correlog(z.slp$x, z.slp$y,
                         z = residuals(error_max_w),
                         na.rm = TRUE,
                         increment = 1,
                         resamp = 1)

plot(cor.lm.res, main="OLS correlation", ylim=c(-3,0.6))
abline(h=0, col="red")
plot(cor.sar1.res, main="SAR correlation",ylim=c(-3,0.6), ylab="", yaxt="n")
abline(h=0, col="red")


# Modelos SAR para cada área ----------------------------------------------

#Voy a correr los modelos paras las áreas transformadas y conservadas por separado y un modelo global 

lm(z ~ slope*cat2,
   data = z.slp)
lm(z ~ slope*cat2,
   data = z.slp) |> 
  summary()

C<-
  z.slp |> 
  filter(cat2=="C")
T<-
  z.slp |> 
  filter(cat2=="T")

lm(z ~ slope,
   data = C) |> 
  summary()

lm(z ~ slope,
   data = T) |> 
  summary()


z.slp |> 
  ggplot(aes(y=z, x=slope, color=cat2))+
  geom_point()+
  #geom_smooth(method="lm", formula=y~x, se=F)+
  geom_abline(aes(intercept=7.02377, slope=1.47603), color="#248f5d") +
  geom_abline(aes(intercept=5.99293, slope=1.76647), color="#f56038") +
  # stat_poly_eq(formula = y~x, 
  #              aes(label = paste(after_stat(eq.label), ..rr.label.., sep = "~~~")), 
  #              parse = TRUE)+
  scale_color_manual(values=c("#248f5d","#f56038"), labels = c("Conserved","Transformed"))+
  labs(x="Terrain slope (median)", y="Species richness", color="Categories") +
  theme_classic()+
  theme(axis.title.y = element_blank())
