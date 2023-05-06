setwd("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/clade_density/clade_density/data")
source("../scripts/findSisterOver.R")
library(phytools)

an<-as.matrix(read.csv("density/cladeDensityStatsAng_100.csv"))[,1]
ge<-as.matrix(read.csv("density/cladeDensityStatsGek_100.csv"))[,1]
ig<-as.matrix(read.csv("density/cladeDensityStatsIgu_100.csv"))[,1]
sc<-as.matrix(read.csv("density/cladeDensityStatsSci_100.csv"))[,1]


ce<-as.matrix(read.csv("density/cladeDensityStatsCet_100.csv"))[,1]
ch<-as.matrix(read.csv("density/cladeDensityStatsChi_100.csv"))[,1]
di<-as.matrix(read.csv("density/cladeDensityStatsDip_100.csv"))[,1]
pr<-as.matrix(read.csv("density/cladeDensityStatsPri_100.csv"))[,1]


pdf("../figs/Figure_3.pdf", width=9, height=6)
vioplot(an,
        ge,
        ig,
        sc,
        ce,
        ch,
        di,
        pr,
        col="red", names=c("Anguimorpha", "Gekkota", "Iguania", "Scincoidea", "Certartiodactyla", "Chiroptera", "Diprotodontia", "Primates"), 
	horizontal=TRUE, xlab="Frequency", las=2)


mean(an)
mean(ge)
mean(sc)
mean(ig)
mean(ce)
mean(ch)
mean(di)
mean(pr)
dev.off()


trAngu<-read.tree("trees/trAngu.tre")
trGekk<-read.tree("trees/trGekk.tree")
trIgua<-read.tree("trees/trIgua.tree")
trScin<-read.tree("trees/trScin.tree")
trCeta<-read.tree("trees/trCeta.tree")
trChir<-read.tree("trees/trChir.tree")
trDipr<-read.tree("trees/trDipr.tree")
trPrim<-read.tree("trees/trPrim.tree")

Oan<-as.matrix(read.csv("overlap/overAngu.csv"))
Oge<-as.matrix(read.csv("overlap/overGekk.csv"))
Oig<-as.matrix(read.csv("overlap/overIgua.csv"))
Osc<-as.matrix(read.csv("overlap/overScin.csv"))

Oce<-as.matrix(read.csv("overlap/overCeta.csv"))
Och<-as.matrix(read.csv("overlap/overChir.csv"))
Odi<-as.matrix(read.csv("overlap/overDipr.csv"))
Opr<-as.matrix(read.csv("overlap/overPrim.csv"))

Oan[is.na(Oan)] <- 0
Oge[is.na(Oge)] <- 0
Oig[is.na(Oig)] <- 0
Osc[is.na(Osc)] <- 0
Oce[is.na(Oce)] <- 0
Och[is.na(Och)] <- 0
Odi[is.na(Odi)] <- 0
Opr[is.na(Opr)] <- 0



densSis <- function(tree, over, dens) {
	final <- numeric()
	for (i in 1:100) {
		x <- findSisterOver(tree[[i]], over)
		names <- data.frame()
		for (j in 1:dim(x)[1]) {
			xx <- strsplit(rownames(x[j, ]), split = "_")[[1]]
			names[j, 1] <- paste(xx[1], xx[2], sep = "_")
			names[j, 2] <- paste(xx[3], xx[4], sep = "_")
		}

		meanDensity <- numeric()
		for (k in 1:dim(x)[1]) {
			meanDensity[k] <- (dens[names[k, 1]] + dens[names[k, 2]])/2
		}
		x <- cbind(x, meanDensity)
		final[i] <- median(subset(x, degrOver > 10, select = "meanDensity")[, 1]) - median(subset(x, degrOver < 10, select = "meanDensity")[, 1])
	}
	final
}

hist(densSis(trAngu, Oan, an), main="Anguimorpha", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trGekk, Oge, ge), main="Gekkota", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trIgua, Oig, ig), main="Iguania", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trScin, Osc, sc), main="Scincoidea", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trCeta, Oce, ce), main="Cetartiodactyla", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trChir, Och, ch), main="Chiroptera", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trDipr, Odi, di), main="Diprotodontia", xlab="Difference in clade density", lty="blank", col="red")
hist(densSis(trPrim, Opr, pr), main="Primates", xlab="Difference in clade density", lty="blank", col="red")
