#Calculation of SES.MNTD
library(picante)
library(iCAMP)

#otu--OTU table with OTUs as rows and samples as columns;
#tree--Phylogenetic tree

otu_match <- match.phylo.comm(tree,otu)
ses_mntd <- NTI.p(otu_match$comm,cophenetic(otu_match$phy),nworker=7)
ses_mntd <- 0 - ses_mntd$NTI
