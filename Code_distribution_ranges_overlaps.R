setwd("~/OneDrive - ufpr.br/projetos/clade_density/data")
library(viridis)

### Code to calculate the frequency distributions of range sizes and range overlap sizes.

an<-as.matrix(read.csv("overlap/overAngu.csv"))
ge<-as.matrix(read.csv("overlap/overGekk.csv"))
ig<-as.matrix(read.csv("overlap/overIgua.csv"))
sc<-as.matrix(read.csv("overlap/overScin.csv"))

ce<-as.matrix(read.csv("overlap/overCeta.csv"))
ch<-as.matrix(read.csv("overlap/overChir.csv"))
di<-as.matrix(read.csv("overlap/overDipr.csv"))
pr<-as.matrix(read.csv("overlap/overPrim.csv"))

an[is.na(an)] <- 0
ge[is.na(ge)] <- 0
ig[is.na(ig)] <- 0
sc[is.na(sc)] <- 0
ce[is.na(ce)] <- 0
ch[is.na(ch)] <- 0
di[is.na(di)] <- 0
pr[is.na(pr)] <- 0

pdf("../figs/Figure_2.pdf", width=9, height=6)
layout(matrix(1:8, ncol=4, byrow=TRUE))
cols<-viridis(3, alpha=0.5)

tx<-"Anguimorpha"
x<-an
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Gekkota"
x<-ge
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Iguania"
x<-ig
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Scincoidea"
x<-sc
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Cetartiodactyla"
x<-ce
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])


tx<-"Chiroptera"
x<-ch
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Diprotodontia"
x<-di
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])

tx<-"Primates"
x<-pr
hist(log(diag(x)), col=cols[1], freq=FALSE, ylim=c(0, 0.18), breaks=seq(0,18, 1), main=tx, lty="blank", xlab="Area log(km)", ylab="Relative frequency")
hist(log(x[upper.tri(x, diag=FALSE)]), freq=FALSE, add=TRUE, col=cols[2], lty="blank")
abline(v=mean(log(diag(x)[diag(x)!=0])),col=cols[1])
abline(v=mean(log(x[upper.tri(x, diag=FALSE)][x[upper.tri(x, diag=FALSE)]!=0])),col=cols[2])


dev.off()
