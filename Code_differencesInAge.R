setwd("~/OneDrive - ufpr.br/projetos/clade_density/data/")
source("../scripts/findSisterOver.R")
library(phytools)


an<-as.matrix(read.csv("overlap/overAngu.csv"))
ge<-as.matrix(read.csv("overlap/overGekk.csv"))
ig<-as.matrix(read.csv("overlap/overIgua.csv"))
sc<-as.matrix(read.csv("overlap/overScin.csv"))

ce<-as.matrix(read.csv("overlap/overCeta.csv"))
ch<-as.matrix(read.csv("overlap/overChir.csv"))
di<-as.matrix(read.csv("overlap/overDipr.csv"))
pr<-as.matrix(read.csv("overlap/overPrim.csv"))

trAngu<-read.tree("trees/trAngu.tre")
trGekk<-read.tree("trees/trGekk.tree")
trIgua<-read.tree("trees/trIgua.tree")
trScin<-read.tree("trees/trScin.tree")

trCeta<-read.tree("trees/trCeta.tree")
trChir<-read.tree("trees/trChir.tree")
trDipr<-read.tree("trees/trDipr.tree")
trPrim<-read.tree("trees/trPrim.tree")


pdf("../figs/Figure_2.pdf", width=9, height=6)
layout(matrix(1:8, ncol=4, byrow=TRUE))

tree<-trAngu
over<-an
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}

hist(difAge, main="Anguimorpha", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trGekk
over<-ge
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Gekkota", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trIgua
over<-ig
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Iguania", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trScin
over<-sc
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Scincoidea", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trCeta
over<-ce
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Cetartiodactyla", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trChir
over<-ch
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Chiroptera", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trDipr
over<-di
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Diprotodontia", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

tree<-trPrim
over<-pr
difAge<-numeric()
for(i in 1:100){
	x<-findSisterOver(tree[[i]], over)
	difAge[i]<-mean(log10(subset(x, degrOver>10, select="mrca")[,1]))-mean(log10(subset(x, degrOver<10, select="mrca")[,1]))
}
hist(difAge, main="Primates", breaks=seq(-0.2534, 0.7, 0.04), xlab="log(Age)", lty="blank", col="red"); abline(v=0, col="blue")
10^mean(difAge); range(difAge)

dev.off()