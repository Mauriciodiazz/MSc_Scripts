


















predict(lm.g, z.slp.sample.Tot, interval = "confidence", level = 0.95)

predict(bes.tot, newdata = z.slp.sample.Tot) |>  # Predicciones del modelo ajustado
  data.frame() |> head()

# Error estándar de las predicciones (aproximación)
se_pred <- sqrt(diag(vcov(bes.tot)))  # Extraer la matriz de varianza-covarianza
alpha <- 0.05  # Nivel de confianza del 95%


z_value <- qnorm(1 - alpha / 2)  # Valor crítico para el 95% (aproximadamente 1.96)
lower_bound <- predictions$trend - z_value * se_pred
upper_bound <- predictions$trend + z_value * se_pred

plot_data <- data.frame(
  x = z.slp.sample.Tot$slope,  # Variable independiente en el eje X
  fit = predictions$trend,       # Predicciones ajustadas
  lower = lower_bound,     # Límite inferior del intervalo de confianza
  upper = upper_bound      # Límite superior del intervalo de confianza
)

ggplot(plot_data, aes(x = x, y = fit)) +
  geom_line(color = "blue", size = 1) +  # Línea ajustada (predicción)
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray") +  # Intervalos de confianza
  labs(title = "Predicción con Intervalos de Confianza",
       x = "Variable Independiente",
       y = "Variable Dependiente") +
  theme_minimal()

