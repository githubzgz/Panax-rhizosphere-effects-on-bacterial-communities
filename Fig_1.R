#The calculation of alpha-diversity and difference rate
library(vegan)

alpha.comput<-function(otu,tree){
  #otu--OTU table with OTUs as rows and samples as columns
  #tree--Phylogenetic tree of bacterial OTUs
  require(vegan)
  N0<-rowSums(otu>0)#OTU richness
  H<-diversity(otu)#Shannon_index
  pd <- pd(otu,tree)$PD
  div<-data.frame(N0,H,pd)
  colnames(div)<-c("Richness","Shannon",'pd')
  return(div)
}

#Difference_rate
alpha_differ <- data.frame(
  difference = ((RS_alpha$Richness - BS_alpha$Richness)/BS_alpha$Richness)*100,
  plant = 'PG',plant_age = 'two')

