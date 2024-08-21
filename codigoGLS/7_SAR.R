
#___________________________________________________________________________#
####                           1_SAr function                 ####
#___________________________________________________________________________#

SARS_fish<- function(coord, tb_sel, vc, name_var) {
  
  #This function operates on coordinates, variable selected, breackboint=vc, and name variable.
  dir0<-getwd()
  dir.create(paste0(gsub(" ","_",name_var))) #create directory with variable name
  setwd(paste0(gsub(" ","_",name_var))) #set directory with variable name
  
  #the function is defined by:
  colnames(coord)<-c("lon","lat")
  coords<-as.data.frame(coord)
  tbtl <- cbind(coords, tb_sel)
  tbtl$absY <- abs(tbtl[,2]) 
  vc <- vc #remeber vC is a unique breakpoint 
  
  
  #split de table using breakpoint in left side rigth side and total table without breakpoint
  tblist <- list()
  tblist[[1]] <- tbtl[which(tbtl$absY < as.numeric(vc)), ] # left table
  tblist[[2]] <- tbtl[which(tbtl$absY >= as.numeric(vc)), ]# Right table
  tblist[[3]] <- tbtl  #All TB
  dir1 <- getwd()
  lf<-list()
  for (n in 1:length(tblist)) {
    
    dir.create(paste0("slopeTB", n))
    setwd(paste0("slopeTB", n))
    
    
    #createcoord
    table <- as.data.frame(tblist[[n]])
    
    sp <- st_as_sf( x = table,coords = c("lon","lat"),
                    crs = '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')#spacialize points
    
    k1 = knn2nb(knearneigh(sp, k = 1)) #create neighbours list of class nb
    
    dist <- unlist(nbdists(k1, sp)) #calculate distances.
    max1 <- max(dist)
    min1 <- min(dist)
    # create a series of neighbour matrices based on different distances (distances are in km)
    
    d1 <- dnearneigh(sp,longlat = F,d1 = 0, d2 = min1)
    d2 <- dnearneigh(sp, longlat = F, d1 = 0,d2 = max1)
    
    ### Create models for each path and weighting scheme ###

    tablelm<-as.data.frame(cbind(table$absY,table[,3]))
    colnames(tablelm)<-c("absY", paste(name_var))
    head(tablelm)
    # We tested three schemes ("W", "C", "S"), two distance (d1, d2) and two models (OLS and SAR)
    scheme <- c("W", "C", "S")
    
    dir2 <- getwd()
    setwd(dir2)
    list0<-list()
    for (i in 1:length(scheme)) {
      # Create a directory by scheme
      dir_scheme <- scheme[i]
      if (!dir.exists(dir_scheme))
        dir.create(dir_scheme)
      setwd(dir_scheme)
      # Create weighting scheme distance matrix
      # Create weighting scheme distance matrix
      set.ZeroPolicyOption(TRUE)
      spatial_weights_d1 <- nb2listw(d1, zero.policy = TRUE, style = scheme[i])
      spatial_weights_d2 <- nb2listw(d2, zero.policy = TRUE, style = scheme[i])
      
      
      schemes <- rep(c(paste(scheme[i])), times = 3)
      mediandef <- rep(paste0(colnames(tablelm)[2], times = 3))
      slopet <- rep(paste0("slope", n), times = 3)
      models <- c("lm_mod", "error_d1", "error_d2")
      results <- data.frame(
        scheme = schemes,
        mediandef = mediandef,
        slopet = slopet,
        model = models,
        Rsquared = NA,
        Rsquared_adj = NA,
        pIntercetp = NA,
        pPendiente = NA,
        Intercept = NA,
        slope = NA,
        Glibert = NA,
        sigma = NA,
        AIC = NA,
        Imoran=NA
      )
      
      # OLS model
      #regression
      lm_mod <-lm(tablelm[, 2] ~ tablelm$absY) # recuerda que para todos los valores de diversificacion ya esta terminando el logaritmo
      summary_reg <- summary(lm_mod)
      #save_coefficients
      results$Rsquared[1] <- summary_reg$r.squared #R-squared
      results$Rsquared_adj[1] <-summary_reg[["adj.r.squared"]] #R-squared adjust
      results$pIntercetp[1] <-summary_reg$coefficients[, "Pr(>|t|)"][1] #p value intercepto
      results$pPendiente[1] <-summary_reg$coefficients[, "Pr(>|t|)"][2] #p value pendiente
      results$Intercept[1] <-summary_reg$coefficients[, "Estimate"][1] #Intercepto
      results$slope[1] <-summary_reg$coefficients[, "Estimate"][2] #pendiente
      results$Glibert[1] <-summary_reg$df[2]                        #grados de libertad
      results$sigma[1] <-summary_reg$sigma                        #sigma_errorstandar
      results$Imoran[1]<-999
      #assing to save
      assign(paste0("reg_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), lm_mod)
      assign(paste0("summary_","slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), summary_reg)
      #assing to save
      capture.output(summary_reg,
                     file = paste0("summary_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]))
      
      
      jpeg(
        filename = paste0("Residuals", n,"_",scheme[i],"_", colnames(tablelm)[2],".jpg"),
        width = 215,
        height = 279,
        units = "mm",
        res = 600
      )
      par(mfrow = c(1, 1))
      plot(tablelm$absY,lm_mod$residuals)
      dev.off()
      
      # SAR error models
      # Minimum distance
      
      error_d1 <- tryCatch({spatialreg::errorsarlm(
        tablelm[, 2] ~ tablelm$absY,
        data = table,
        listw = spatial_weights_d1,
        tol = 1e-20,
        zero.policy = T
      )
        # Code that might generate an error goes here
      }, error = function(e) {
        
        error_d1 <-spatialreg::errorsarlm(
          tablelm[, 2] ~ tablelm$absY,
          data = table,
          listw = spatial_weights_d1,
          tol = 1e-20,
          zero.policy = T)
        
        # Code to handle the error goes here, if desired
      })
      
      
      error_d1_s <- summary(error_d1, Nagelkerke = TRUE)
      R2_error_d1 <- error_d1_s$NK
      var<- tablelm[, 2]
      Id1<-moran(var,  spatial_weights_d1, length(k1), Szero(spatial_weights_d1))
      
      results$Rsquared[2] <- error_d1_s$NK #R-squared
      results$Rsquared_adj[2] <- error_d1_s$NK #R-squared adjust
      results$pIntercetp[2] <-error_d1_s$Coef[, "Pr(>|z|)"][1]                  #p value intercepto
      results$pPendiente[2] <-error_d1_s$Coef[, "Pr(>|z|)"][2]                  #p value pendiente
      results$Intercept[2] <- error_d1_s$Coef[, "Estimate"][1]                  #Intercepto
      results$slope[2] <- error_d1_s$Coef[, "Estimate"][2]                      #pendiente
      results$Glibert[2] <-length(error_d1_s$y) - error_d1_s$parameters         #grados de libertad
      results$sigma[2] <- error_d1_s$lambda.se                                  #sigma_lamnda
      results$Imoran[2]<-Id1$I
      
      #asigno para guardarla
      assign(paste0("modelE1", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), error_d1)
      assign(paste0("summaryE1_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), error_d1_s)
      #guardo texto salida
      capture.output(error_d1_s,
                     file = paste0("summary_E1_","slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]))
      
      # Maximum distance
      
      error_d2 <- tryCatch({spatialreg::errorsarlm(
        tablelm[, 2] ~ tablelm$absY,
        data = table,
        listw = spatial_weights_d2,
        tol = 1e-12,
        zero.policy = T
      )
        # Code that might generate an error goes here
      }, error = function(e) {
        
        error_d2 <-spatialreg::errorsarlm(
          tablelm[, 2] ~ tablelm$absY,
          data = table,
          listw = spatial_weights_d2,
          tol = 1e-20,
          zero.policy = T)
        
        # Code to handle the error goes here, if desired
      })
      
      
      error_d2_s <- summary(error_d2, Nagelkerke = TRUE)
      R2_error_d2 <-error_d2_s$NK#
      var<- tablelm[, 2]
      Id2<- moran(var,  spatial_weights_d2, length(k1), Szero(spatial_weights_d2))
      
      results$Rsquared[3] <- error_d2_s$NK #R-squared
      results$Rsquared_adj[3] <- error_d2_s$NK #R-squared adjust
      results$pIntercetp[3] <-error_d2_s$Coef[, "Pr(>|z|)"][1] #p value intercepto
      results$pPendiente[3] <-error_d2_s$Coef[, "Pr(>|z|)"][2] #p value pendiente
      results$Intercept[3] <-    error_d2_s$Coef[, "Estimate"][1] #Intercepto
      results$slope[3] <- error_d2_s$Coef[, "Estimate"][2] #pendiente
      results$Glibert[3] <- length(error_d2_s$y) - error_d2_s$parameters                      #grados de libertad
      results$sigma[3] <- error_d2_s$lambda.se    #sigma_lamnda
      results$Imoran[3]<-Id2$I
      #asigno para guardarla
      assign(paste0("modelE2_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), error_d2)
      assign(paste0("summaryE2_","slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]), error_d2_s)
      #guardo texto salida
      capture.output(error_d2_s,
                     file = paste0("summary_E2_","slope_", n,"_",scheme[i],"_", colnames(tablelm)[2]),error_d2_s)
      
      
      # Save R2 (pseudo R2 for errorsar) and AIC
      results$AIC <- AIC(lm_mod, error_d1, error_d2)
      
      write.csv(results, file = paste0("slope_", n,"_",scheme[i],"_", colnames(tablelm)[2],".csv"),
                row.names = F)
      # Make correlograms of residual autocorrelation
      cor.ols1.res <-
        correlog(tablelm[, 2],tablelm$absY,
                 z = residuals(lm_mod),
                 na.rm = TRUE,
                 increment = 1,
                 resamp = 1
        )
      cor.sar1.res <-
        correlog(tablelm[, 2],tablelm$absY,
                 z = residuals(error_d1),
                 na.rm = TRUE,
                 increment = 1,
                 resamp = 1
        )
      cor.sar2.res <-
        correlog(tablelm[, 2],tablelm$absY,
                 z = residuals(error_d2),
                 na.rm = TRUE,
                 increment = 1,
                 resamp = 1
        )
      # Save correlograms in a jpeg file
      jpeg(
        filename = paste0("slope_", n,"_",scheme[i],"_", colnames(tablelm)[2],".jpg"),
        width = 215,
        height = 279,
        units = "mm",
        res = 600
      )
      par(mfrow = c(3, 1))
      plot(
        cor.ols1.res,
        xlab = "Distance (Km)",
        ylab = "Moran's I",
        ylim = c(-1, 1),
        type = "l",
        lwd = 2,
        main = paste(models[1]),
        cex.main = 2,
        cex.lab = 1.8,
        cex.axis = 1.5
      )
      abline(h = 0, lty = 5)
      plot(
        cor.sar1.res,
        xlab = "Distance (Km)",
        ylab = "Moran's I",
        ylim = c(-1, 1),
        type = "l",
        lwd = 2,
        main = paste(models[2]),
        cex.main = 2,
        cex.lab = 1.8,
        cex.axis = 1.5
      )
      abline(h = 0, lty = 5)
      plot(
        cor.sar2.res,
        xlab = "Distance (Km)",
        ylab = "Moran's I",
        ylim = c(-1, 1),
        type = "l",
        lwd = 2,
        main = paste(models[3]),
        cex.main = 2,
        cex.lab = 1.8,
        cex.axis = 1.5
      )
      abline(h = 0, lty = 5)
      dev.off()
      # Save path results
      List_data2<-list(
        sarE2= get(paste0("modelE2_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2])),
        sarE1= get(paste0("modelE1", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2])),
        ols1=get( paste0("reg_", "slope_", n,"_",scheme[i],"_", colnames(tablelm)[2])),
        corols=cor.ols1.res ,
        corsar1= cor.sar1.res, 
        corsar2=  cor.sar2.res) 
      
      assign(paste0(scheme[i]),List_data2)
      setwd(dir2)
    }
    Scheme<-list(W=W,C=C,S=S)
    assign(paste0("slope",n),Scheme)
    setwd(dir1)
    
  }
  lf<-list(slope1=slope1,slope2=slope2,slope3=slope3)
  
  setwd(dir0)
  return(lf)
}


#___________________________________________________________________________#
####                           2_SAr for ALLDR           ####
#___________________________________________________________________________#

load(here::here("AssamblajeLevel/OLS/OLSGAM.RData"))
setwd(here::here("AssamblajeLevel"))
dir.create("SAR")
setwd("SAR")
dir.create("DR")
setwd("DR")

coordTot<-comparetbTot[,c(1,2)]
dim(Speciations)

directo<-getwd()
for (x in 1:dim(Speciations)[2]) {
  setwd(directo) 
  print(x) 
  coord<-coordTot
  tb_sel<-Speciations[,x]
  name_var<-colnames(Speciations)[x]
  vc<-DRvc[x]
  olssar<-SARS_fish(coord, tb_sel,vc, name_var)
  assign(paste0("olssar",gsub(" ","_",name_var)),olssar)
  save(olssar,file=paste0("olssar",gsub(" ","_",name_var),".RData"))
  
}


####Tb_resume ####



folders <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)

csvfiles<-c()
for (u in 1:length(folders)) {
  files<-list.files(path = paste0(folders[u]), '.csv$', full.names = TRUE)
  csvfiles<-c(csvfiles,files) 
}

tb_slopes<- do.call("rbind", lapply(csvfiles[grep("slope", csvfiles )], read.csv))


write.csv(tb_slopes,file="FinalReportDR.csv")





#___________________________________________________________________________#
####                           1_SAr for SubPAM           ####
#___________________________________________________________________________#
setwd("../")
dir.create("BAMDRCLADS")
setwd("BAMDRCLADS")

directo<-getwd()
coordTot<-comparetbSub[,c(1,2)]
colnames(coordTot)
save.image("OLSSAR_BAMCLDR.RData")
for (x in 1:dim(Speciations)[2]) {
  setwd(directo) 
  print(x) 
  coord<-coordTot
  tb_sel<-Speciations[,x]
  name_var<-colnames(Speciations)[x]
  vc<-vcBAM[x]
  olssar<-SARS_fish(coord, tb_sel,vc2, name_var)
  assign(paste0("olssar",gsub(" ","_",name_var)),olssar)
  save(olssar,file=paste0("olssar",gsub(" ","_",name_var),".RData"))
  
}

save.image(resultSAR.RData)


##Tbresumen##


folders <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)

csvfiles<-c()
for (u in 1:length(folders)) {
  files<-list.files(path = paste0(folders[u]), '.csv$', full.names = TRUE)
  csvfiles<-c(csvfiles,files) 
}

tb_slopes<- do.call("rbind", lapply(csvfiles[grep("slope", csvfiles )], read.csv))


write.csv(tb_slopes,file="FinalReportBA.csv")







