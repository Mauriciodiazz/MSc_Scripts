
data(iris)
attach(iris)
length(Sepal.Length)
plot(Sepal.Length, Petal.Length)
Miris <- lm(Petal.Length~Sepal.Length, data=iris)
Miris
summary(Miris)
anova(Miris)
par(mfrow=c(2,2))
plot(Miris) # aunque el modelo es significativo, el gráfico de residuales vs ajustados muestra que el modelo no es muy ajustado. los residuales no muestran homocedasticidad en los valores grandes. ES NECESARIO BUSCAR OTRO MODELO


plot(Sepal.Length, Petal.Length)
abline(Miris, col="blue")
predict_interval <- predict(Miris, interval="confidence", level=0.95)
lines(Sepal.Length, predict_interval[, "lwr"], col="red", lty=2)
lines(Sepal.Length, predict_interval[, "upr"], col="red", lty=2)

# para este caso, los datos provienen de 3 especies, lo que puede estar causando un desajuste del modelo

Miris2 <- lm(Petal.Length~Sepal.Length*Species) #incluyendo la especie
par(mfrow=c(2,2))
plot(Miris2)
Miris2
summary(Miris2)

plot(Sepal.Length, Petal.Length)
abline(h=mean(Petal.Length), col="blue") # hipótesis nula

species_names <- unique(iris$Species)  #hipótesis nulas de cada spp
for (species in species_names) {
  mean_length <- mean(iris$Petal.Length[iris$Species == species])
  abline(h=mean_length, col="red")
}

#construir vectores para mantener unicamente la información que necesitamos para el cálculo de interceptos y pendientes para cada especie, con el orden en que aparecen los coeficientes en el modelo (en este caso Miris2)

ISe <- c(1,0,0,0,0,0)
IVe <- c(1,0,1,0,0,0)
IVir <- c(1,0,0,1,0,0)

Pse <- c(0,1,0,0,0,0)
Pve <- c(0,1,0,0,1,0)
PVir <- c(0,1,0,0,0,1)
Pgral<- (Pse+Pve+PVir)/3 # para calcular la pendiente general, que es el promedio de las 3 pendientes
Pve-PVir


library(gmodels)

estimable(Miris2, rbind(Pse, Pve, PVir))
estimable(Miris2, rbind(se=Pse, ve=Pve, vir=PVir, se_ve=Pse-Pve, se_vir=Pse-PVir, ve_vir=Pve-PVir))

plot(Sepal.Length, Petal.Length, las=1, type = "n")
points(Sepal.Length[Species =="setosa"], Petal.Length[Species =="setosa"], col="red")
points(Sepal.Length[Species =="versicolor"], Petal.Length[Species =="versicolor"], col="limegreen")
points(Sepal.Length[Species =="virginica"], Petal.Length[Species =="virginica"], col="black")

qt(0.975, 144, lower.tail = F) #para encontrar el valor que multiplica la predicción para los intervalos de confianza

SEQ <- seq(min(Sepal.Length[Species =="setosa"]), max(Sepal.Length[Species =="setosa"]), length=30) #establecer valores entre el mínimio y el máximo de las observaciones
PRED <- predict(Miris2, data.frame(Sepal.Length = SEQ, Species ="setosa"), se = T) # "se" debe ponerse para incluir el error estandar
lines(SEQ, PRED$fit, col="darkred", lwd=2)# línea de regresión
lines(SEQ, PRED$fit-PRED$se.fit*1.976575, col="darkred", lwd=2, lty=2) #intervalo de confianza inferior
lines(SEQ, PRED$fit+PRED$se.fit*1.976575, col="darkred", lwd=2, lty=2) #intervalo de confianza superior



SEQ <- seq(min(Sepal.Length[Species =="versicolor"]), max(Sepal.Length[Species =="versicolor"]), length=30)
PRED <- predict(Miris2, data.frame(Sepal.Length = SEQ, Species ="versicolor"), se = T)
lines(SEQ, PRED$fit, col="limegreen", lwd=2)
lines(SEQ, PRED$fit-PRED$se.fit*1.976575, col="limegreen", lwd=2, lty=2)
lines(SEQ, PRED$fit+PRED$se.fit*1.976575, col="limegreen", lwd=2, lty=2)



SEQ <- seq(min(Sepal.Length[Species =="virginica"]), max(Sepal.Length[Species =="virginica"]), length=30)
PRED <- predict(Miris2, data.frame(Sepal.Length = SEQ, Species ="virginica"), se = T)
lines(SEQ, PRED$fit, col="black", lwd=2)
lines(SEQ, PRED$fit-PRED$se.fit*1.976575, col="black", lwd=2, lty=2)
lines(SEQ, PRED$fit+PRED$se.fit*1.976575, col="black", lwd=2, lty=2)




abline(h=mean(Petal.Length[Species=="virginica"]), col="black") #h0
abline(h=mean(Petal.Length[Species=="setosa"]), col="red")#h0
abline(h=mean(Petal.Length[Species=="versicolor"]), col="limegreen")#h0


summary(Miris2)
