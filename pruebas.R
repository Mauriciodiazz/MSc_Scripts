


sp_g <- read_sf("F:/Maestria_DD/Shapes_MSc_DD/borrar/prueba_SAR.shp")
  # st_as_sf( x = z.slp, 
  #                 coords = c("x","y"),
  #                 crs = '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')
write_sf(sp_g,"F:/Maestria_DD/Shapes_MSc_DD/borrar/prueba_SAR2.shp")

k1_g = knn2nb(knearneigh(sp_g, k = 1)) #establece cuáles son los vecinos según el k. El warning es porque hay puntos que no tienen los 5 vecinos (pe, las esquinas)

#dist_g <- 
  nbdists(k1_g, sp_g, longlat = T) #|> unlist()

min1_g <- min(dist_g)

d_min_g <- dnearneigh(sp_g, longlat = F, d1 = 0, d2 = min1_g)

sw_min_g <- nb2listw(d_min_g, zero.policy = TRUE, style = "W") 

# SAR
error_min_g <- spatialreg::errorsarlm(
  z ~ slope,
  data = z.slp,
  listw = sw_min_g,
  tol = 1e-12, 
  zero.policy = TRUE)

summary(error_min_g, Nagelkerke=TRUE)

plot(z.slp$slope, z.slp$z)



lets.correl(residuals(lm.g), 
            lets.distmat(as.matrix(z.slp[,1:2])), 12,
            equidistant = FALSE, 
            plot = TRUE)
lets.correl(residuals(error_min_g), 
            lets.distmat(as.matrix(z.slp[,1:2])), 12,
            equidistant = FALSE, 
            plot = TRUE)


lets.correl(residuals(lm.c), 
            lets.distmat(as.matrix(C[,1:2])), 12,
            equidistant = FALSE, 
            plot = TRUE)
lets.correl(residuals(error_min_c), 
            lets.distmat(as.matrix(C[,1:2])), 12,
            equidistant = FALSE, 
            plot = TRUE)


correlog(z.slp$x, z.slp$y,
         z = residuals(error_min_g),
         na.rm = TRUE,
         increment = 1,
         resamp = 1)

