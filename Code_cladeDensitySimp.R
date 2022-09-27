rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/clade_densi")

library(phytools)

trCe<-read.tree("data/trCeta.tree")
trAn<-read.nexus("data/anguimorpha_tree.nex")

datCe<-as.matrix(read.csv("data/Ceta_overlap-triang-matrix.csv", row.names=1,
                          header=T))
datAn<-as.matrix(read.csv("data/ANGUI_overlap.csv", row.names=1,
                          header=T))

datCe<-apply(datCe, c(1,2), as.numeric)
datAn<-apply(datAn, c(1,2), as.numeric)

datCe[is.na(datCe)] <- 0
datAn[is.na(datAn)] <- 0

cladeDensity<-function(mat, tree, log=TRUE){
  library(ape)
  library(geiger)
  source("https://raw.githubusercontent.com/mgharvey/ES-sim/master/R/essim.R")
  
  mat<-mat[order(rownames(mat)),]
  mat<-mat[,order(colnames(mat))]
  if(log==TRUE){
    mat<-mat+1
    mat<-log10(mat)
  }
  x<-mat #/diag(mat)

  if(missing(tree)){
    diag(x)<-0
    cladeDen<-rowSums(x)
    res<-data.frame(Species=names(cladeDen), 
                    Clade_density=cladeDen)
    colnames(res)<-c("Species", "Clade_density")
    rownames(res)<-NULL
    return(res)
  } else {
    
    es<-compute_es(tree)
    
    tr<-treedata(tree, mat, warnings=F)$phy
    
    phyl<-vcv(tr)
    phyl<-phyl[order(rownames(phyl)),]
    phyl<-phyl[,order(colnames(phyl))]
    
    check1<-setdiff(rownames(phyl), rownames(mat))
    check2<-setdiff(rownames(mat), rownames(phyl))
    check<-list(sort(check1), sort(check2))
    names(check)<-cbind("Taxa not found in data",
                        "Taxa not found in tree")
    if (length(check1)==0 && length(check2)==0) {
      x<-x*phyl
      diag(x)<-0
      cladeDen<-rowSums(x)
      
      es<-es[names(es) %in% names(cladeDen)]
      es<-es[order(names(es))]
      
      resEssim<-essim(phy=tr, trait=cladeDen, nsim=100, es=es)
      
      res<-data.frame(Species=names(cladeDen), 
                      Clade_density=cladeDen, es=es)
      colnames(res)<-c("Species", "Clade_density", "es")
      rownames(res)<-NULL
      final<-list()
      final$CD<-res
      final$essim<-resEssim
      return(final)
    } else {
      return(check)
    }
  }
}

tr<-trCe
dat<-datCe
claDenLog<-lapply(tr, cladeDensity, mat=dat, log=T)

CD<-matrix(ncol=length(tr), nrow=nrow(claDenLog[[1]]$CD))
for(i in 1:length(tr)){
  CD[,i]<-claDenLog[[i]]$CD[,2]
}
rownames(CD)<-claDenLog[[1]]$CD[,1]

Mean<-Ci_inf<-Ci_sup<-numeric()
for(i in 1:nrow(CD)){
  Mean[i]<-mean(CD[i,])
  Ci_inf[i]<-quantile(CD[i,], probs=0.025)
  Ci_sup[i]<-quantile(CD[i,], probs=0.975)
}
CDstats<-cbind(Mean, Ci_inf, Ci_sup)
rownames(CDstats)<-rownames(CD)
write.csv(CDstats, "data/cladeDensityStatsCet.csv")



tr<-trAn
dat<-datAn
claDenLog<-lapply(tr, cladeDensity, mat=dat, log=T)

CD<-matrix(ncol=length(tr), nrow=nrow(claDenLog[[1]]$CD))
for(i in 1:length(tr)){
  CD[,i]<-claDenLog[[i]]$CD[,2]
}
rownames(CD)<-claDenLog[[1]]$CD[,1]

Mean<-Ci_inf<-Ci_sup<-numeric()
for(i in 1:nrow(CD)){
  Mean[i]<-mean(CD[i,])
  Ci_inf[i]<-quantile(CD[i,], probs=0.025)
  Ci_sup[i]<-quantile(CD[i,], probs=0.975)
}
CDstats<-cbind(Mean, Ci_inf, Ci_sup)
rownames(CDstats)<-rownames(CD)
write.csv(CDstats, "data/cladeDensityStatsAng.csv")

