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

mod.null= glm(z ~ 1,
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



