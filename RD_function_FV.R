### Function to calculate the RangeDiversity parameters (based on Arita et al. 2012. GEB)

### DATA: 
### Delta = a presence-absence matrix (PAM); ***species are rows and sites are columns*** (without coordinates or any other information)
### *Note that most PAMs nowadays come as transpose; with species as columns and sites as rows. 
### *If so, you'd need to transpose your PAM before using this function. 

RangeDiversity <- function(Delta){
  
  DeltaT <- t(Delta)
  
  #Calculate parameters of range and richness
  
  Species<- nrow(Delta)
  Quadrats<- ncol(Delta)
  One_S<-as.matrix(seq(1, 1, length=Species))
  One_N<-as.matrix(seq(1, 1, length=Quadrats))
  
  SpeciesRichness<- as.matrix(colSums(Delta))
  RangeSize<- as.matrix(rowSums(Delta))
  
  Fill<- sum(RangeSize)
  Beta<- Species*Quadrats/Fill
  
  RangeMean<- mean(RangeSize)
  RangeMin<- min(RangeSize)
  RangeMax<- max(RangeSize)
  D_Volume<- Delta%*%SpeciesRichness
  RangeRichness<- D_Volume/RangeSize
  SpeciesCovariance<-(Delta%*%DeltaT-RangeSize%*%t(RangeSize)/Quadrats)/Quadrats
  RangeVariance<- RangeSize/Quadrats*(1-RangeSize/Quadrats)
  VarianceOfRanges<-t(RangeSize)%*%RangeSize/Species-(t(One_S)%*%RangeSize/Species)^2
  
  RichnessMean<- mean(SpeciesRichness)
  RichnessMin<- min(SpeciesRichness)
  RichnessMax<- max(SpeciesRichness)
  R_Volume<- DeltaT%*%RangeSize
  SiteRange<- R_Volume/SpeciesRichness
  SitesCovariance<- (DeltaT%*%Delta-SpeciesRichness%*%t(SpeciesRichness)/Species)/Species
  RichnessVariance<- SpeciesRichness/Species*(1-SpeciesRichness/Species)
  VarianceOfRichness<-t(SpeciesRichness)%*%SpeciesRichness/Quadrats-(t(One_N)%*%SpeciesRichness/Quadrats)^2
  
  CovarianceSpecies<- RangeSize*(RangeRichness-RichnessMean)/Quadrats/Species
  CovarianceSites<- SpeciesRichness*(SiteRange-RangeMean)/Species/Quadrats
  
  V_species<- VarianceOfRichness/sum(RangeVariance)
  V_sites<- VarianceOfRanges/sum(RichnessVariance)
  U_species<- sum(abs(RangeSize-RangeMean)/Quadrats)/Species
  U_sites<- sum(abs(SpeciesRichness-RichnessMean)/Species)/Quadrats
  
  list(RD.single.values=data.frame(Species, Quadrats, RichnessMean, RangeMean,Fill, Beta, V_species, V_sites),RD.species.values=data.frame(RangeSize,RangeRichness,RangeVariance,CovarianceSpecies),RD.sites.values=data.frame(SpeciesRichness,SiteRange,RichnessVariance,CovarianceSites),SitesCovariance=SitesCovariance,SpeciesCovariance=SpeciesCovariance)  
    
}
