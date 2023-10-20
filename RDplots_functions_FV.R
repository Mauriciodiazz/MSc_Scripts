#### by Species ####

### DATA: 
### pamRD is the output from applying the "RangeDiversity" function (file: RD_function.R)

species.RDplot<- function(pamRD) {
  
  Species <- pamRD$RD.single.values$Species
  Quadrats <- pamRD$RD.single.values$Quadrats
  RangeSize <- pamRD$RD.species.values$RangeSize
  RangeRichness <- pamRD$RD.species.values$RangeRichness
  RichnessMin <- min(pamRD$RD.sites.values$SpeciesRichness)
  RichnessMax <- max(pamRD$RD.sites.values$SpeciesRichness)
  RichnessMean <- mean(pamRD$RD.sites.values$SpeciesRichness)
  
  def.par <-par(no.readonly = TRUE)
  xhist<- hist(RangeRichness/Species, breaks = seq(0, 1, .05), plot = FALSE)
  yhist<- hist(RangeSize/Quadrats, breaks = seq(0, 1, .05), plot = FALSE)
  topx<-max(xhist$counts)
  topy<-max(yhist$counts)
  nf<- layout(matrix(c(2, 0, 1, 3), 2, 2, T), c(6, 1), c(1, 6), TRUE)
  layout.show(nf)
  
  #RD plot
  par(mar= c(5, 5, 1, 1))
  plot(RangeRichness/Species, RangeSize/Quadrats, xlim = c(0,1), ylim = c(0, 1), xlab = "Mean proportional range richness", ylab = "Proportional range size",pch=21,bg="black",cex=1.5)
  
  
  #Isocovariance lines
  x<- seq(RichnessMean/Species, 1, length = 100)
  for(i in c(.01, .05, .1)) {
    lines(x, i/(x- RichnessMean/Species), lwd = 2, col = "pink")
  }
  y<- seq(0, RichnessMean/Species, length = 100)
  for(j in c(-.01, -.05, -.1)) {
    lines(y, j/(y-RichnessMean/Species), lwd = 2, col = "pink")
  }
  
  
  x<- seq(RichnessMin/Species, RichnessMean/Species, length = 100)
  lines(x, (RichnessMax-RichnessMean)/Species/(RichnessMax/Species-x), lwd = 2)
  y<- seq(RichnessMean/Species, RichnessMax/Species, length = 100)
  lines(y, (RichnessMean-RichnessMin)/Species/(y-RichnessMin/Species), lwd = 2)
  segments(RichnessMean/Species, 0, RichnessMean/Species, 1, lty = 3, lwd = 1.5)
  
  #Top histogram
  par(mar=c(0, 5, 1, 1))
  barplot(xhist$counts, axes= FALSE, ylim = c(0, topx), space =0)
  
  #Side histogram
  par(mar=c(5, 0, 1, 1))
  barplot(yhist$counts, axes = FALSE, space = 0, horiz = TRUE)
  
  par(def.par)
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
  
  
  #ISOCOVARIANCE LINES
  x<- seq(RangeMean/Quadrats, 1, length = 100)
  for(i in c(.01, .05, .1)) {
    lines(x, i/(x - RangeMean/Quadrats), lwd = 2, col = "pink")
  }
  y<- seq(0, RangeMean/Quadrats, length = 100)
  for(j in c(-.01, -.05, -.1)) {
    lines(y, j/(y- RangeMean/Quadrats), lwd = 2, col = "pink")
  }

  
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
