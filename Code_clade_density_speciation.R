rm(list = ls())

library(phytools)
library(ape)
library(geiger)
source("https://raw.githubusercontent.com/mgharvey/ES-sim/master/R/essim.R")

### Code to estimate the frequency distribution of the slopes of relationships between interspecific variation in clade density and speciation rate (λDR) across different tested topologies in each taxon.

datAn <- as.matrix(read.csv("ANGUI_overlap.csv", row.names = 1, header = T))
datCe <- as.matrix(read.csv("Ceta_overlap-triang-matrix.csv", row.names = 1,
                          header = T))
datCh <- as.matrix(read.csv("CHIRO_overlap.csv", row.names = 1, header = T))
datDi <- as.matrix(read.csv("Dipro_overlap-triang-matrix.csv", row.names = 1,
                          header = T))
datPr <- as.matrix(read.csv("Primate_overlap-triang-matrix.csv", row.names = 1,
                          header = T))
datSc <- as.matrix(read.csv("SCINC_overlap.csv", row.names = 1, header = T))
datGe <- as.matrix(read.csv("GEKK_overlap.csv", row.names = 1, header = T))
datIg <- as.matrix(read.csv("overIgua.csv", row.names = 1, header = T))

for(i in 1:nrow(datCe)) {
  x <- strsplit(rownames(datCe)[i], " ")[[1]]
  rownames(datCe)[i] <- paste0(x[1], "_", x[2])
  y <- strsplit(colnames(datCe)[i], "[.]")[[1]]
  colnames(datCe)[i] <- paste0(y[1], "_", y[2])
}

for(i in 1:nrow(datCh)) {
  x <- strsplit(rownames(datCh)[i], " ")[[1]]
  rownames(datCh)[i] <- paste0(x[1], "_", x[2])
  y <- strsplit(colnames(datCh)[i], "[.]")[[1]]
  colnames(datCh)[i] <- paste0(y[1], "_", y[2])
}

for(i in 1:nrow(datDi)) {
  x <- strsplit(rownames(datDi)[i], " ")[[1]]
  rownames(datDi)[i] <- paste0(x[1], "_", x[2])
  y <- strsplit(colnames(datDi)[i], "[.]")[[1]]
  colnames(datDi)[i] <- paste0(y[1], "_", y[2])
}

for(i in 1:nrow(datPr)) {
  x <- strsplit(rownames(datPr)[i], " ")[[1]]
  rownames(datPr)[i] <- paste0(x[1], "_", x[2])
  y <- strsplit(colnames(datPr)[i], "[.]")[[1]]
  colnames(datPr)[i] <- paste0(y[1], "_", y[2])
}

for(i in 1:nrow(datGe)) {
  x <- strsplit(rownames(datGe)[i], " ")[[1]]
  rownames(datGe)[i] <- paste0(x[1], "_", x[2])
  y <- strsplit(colnames(datGe)[i], "[.]")[[1]]
  colnames(datGe)[i] <- ifelse(length(y)  =  =  2,paste0(y[1], "_", y[2]),
                             paste0(y[1], "_", y[2], "-", y[3]))
}

datAn <- apply(datAn, c(1,2), as.numeric)
datCe <- apply(datCe, c(1,2), as.numeric)
datCh <- apply(datCh, c(1,2), as.numeric)
datDi <- apply(datDi, c(1,2), as.numeric)
datPr <- apply(datPr, c(1,2), as.numeric)
datSc <- apply(datSc, c(1,2), as.numeric)
datGe <- apply(datGe, c(1,2), as.numeric)
datIg <- apply(datIg, c(1,2), as.numeric)

datAn[is.na(datAn)] <- 0
datCe[is.na(datCe)] <- 0
datCh[is.na(datCh)] <- 0
datDi[is.na(datDi)] <- 0
datPr[is.na(datPr)] <- 0
datSc[is.na(datSc)] <- 0
datGe[is.na(datGe)] <- 0
datIg[is.na(datIg)] <- 0

tAn <- t(datAn)
tCe <- t(datCe)
tCh <- t(datCh)
tDi <- t(datDi)
tPr <- t(datPr)
tSc <- t(datSc)
tGe <- t(datGe)
tIg <- t(datIg)

overAn <- tAn + datAn
overCe <- tCe + datCe
overCh <- tCh + datCh
overDi <- tDi + datDi
overPr <- tPr + datPr
overSc <- tSc + datSc
overGe <- tGe + datGe
overIg <- tIg + datIg

diag(overAn) <- 0
diag(overCe) <- 0
diag(overCh) <- 0
diag(overDi) <- 0
diag(overPr) <- 0
diag(overSc) <- 0
diag(overGe) <- 0
diag(overIg) <- 0

## Trees from VertLife
trMa <- read.nexus("mammalia_node_dated_27JUL2020.nex")[1:100]
trSq <- read.nexus("squamata27JUL20.nex")[1:100]

cladeDensity <- function(mat, tree, log  =  TRUE) {
  
  mat <- mat[order(rownames(mat)), ]
  mat <- mat[, order(colnames(mat))]
  if(log  =  =  TRUE) {
    mat <- mat + 1
    mat <- log10(mat)
  }
  x <- mat

  if(missing(tree)) {
    diag(x) <- 0
    cladeDen <- rowSums(x)
    res <- data.frame(Species  =  names(cladeDen), Clade_density  =  cladeDen)
    colnames(res) <- c("Species", "Clade_density")
    rownames(res) <- NULL
    return(res)
  } else {
    
    es <- compute_es(tree)
    
    tr <- treedata(tree, mat, warnings  =  F)$phy
    
    phyl <- vcv(tr)
    phyl <- phyl[order(rownames(phyl)), ]
    phyl <- phyl[, order(colnames(phyl))]
    
    check1 <- setdiff(rownames(phyl), rownames(mat))
    check2 <- setdiff(rownames(mat), rownames(phyl))
    check <- list(sort(check1), sort(check2))
    names(check) <- cbind("Taxa not found in data", "Taxa not found in tree")
    if (length(check1)  =  =  0 && length(check2)  =  =  0) {
      x <- x * phyl
      diag(x) <- 0
      cladeDen <- rowSums(x)
      
      es <- es[names(es) %in% names(cladeDen)]
      es <- es[order(names(es))]
      
      resEssim <- essim(phy  =  tree, trait  =  cladeDen, nsim  =  100, es  =  es)
      
      res <- data.frame(Species  =  names(cladeDen), 
                        Clade_density  =  cladeDen, es  =  es)
      colnames(res) <- c("Species", "Clade_density", "es")
      rownames(res) <- NULL
      final <- list()
      final$CD <- res
      final$essim <- resEssim
      return(final)
    } else {
      return(check)
    }
  }
}

CDAn <- matrix(ncol  =  length(trSq), nrow  =  nrow(overAn))
essimAn  <-  matrix(ncol  =  2, nrow  =  length(trSq))
for(i in 1:length(trSq)) {
  claDenLogAn <- cladeDensity(mat = overAn, tree = trSq[[i]], log = T)
  CDAn[,i] <- claDenLogAn$CD[,2]
  essimAn[i,] <- claDenLogAn$essim[1:2]
}
rownames(CDAn) <- claDenLogAn$CD[,1]
colnames(essimAn) <- c("rho", "P")

MeanAn <- Ci_infAn <- Ci_supAn <- numeric()
for(i in 1:nrow(CDAn)) {
  MeanAn[i] <- mean(CDAn[i,])
  Ci_infAn[i] <- quantile(CDAn[i,], probs = 0.025)
  Ci_supAn[i] <- quantile(CDAn[i,], probs = 0.975)
}
CDstatsAn <- cbind(MeanAn, Ci_infAn, Ci_supAn)
rownames(CDstatsAn) <- rownames(CDAn)
colnames(CDstatsAn) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsAn, "cladeDensityStatsAng_100.csv")
write.csv(essimAn, "essimAng_100.csv")


CDCe <- matrix(ncol = length(trMa), nrow = nrow(overCe))
essimCe <- matrix(ncol = 2, nrow = length(trMa))
for(i in 1:length(trMa)) {
  claDenLogCe <- cladeDensity(mat = overCe, tree = trMa[[i]], log = T)
  CDCe[,i] <- claDenLogCe$CD[,2]
  essimCe[i,] <- claDenLogCe$essim[1:2]
}
rownames(CDCe) <- claDenLogCe$CD[,1]
colnames(essimCe) <- c("rho", "P")

MeanCe <- Ci_infCe <- Ci_supCe <- numeric()
for(i in 1:nrow(CDCe)) {
  MeanCe[i] <- mean(CDCe[i,])
  Ci_infCe[i] <- quantile(CDCe[i,], probs = 0.025)
  Ci_supCe[i] <- quantile(CDCe[i,], probs = 0.975)
}
CDstatsCe <- cbind(MeanCe, Ci_infCe, Ci_supCe)
rownames(CDstatsCe) <- rownames(CDCe)
colnames(CDstatsCe) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsCe, "cladeDensityStatsCet_100.csv")
write.csv(essimCe, "essimCet_100.csv")


CDCh <- matrix(ncol = length(trMa), nrow = nrow(overCh))
essimCh <- matrix(ncol = 2, nrow = length(trMa))
for(i in 1:length(trMa)) {
  claDenLogCh <- cladeDensity(mat = overCh, tree = trMa[[i]], log = T)
  CDCh[,i] <- claDenLogCh$CD[,2]
  essimCh[i,] <- claDenLogCh$essim[1:2]
}
rownames(CDCh) <- claDenLogCh$CD[,1]
colnames(essimCh) <- c("rho", "P")

MeanCh <- Ci_infCh <- Ci_supCh <- numeric()
for(i in 1:nrow(CDCh)) {
  MeanCh[i] <- mean(CDCh[i,])
  Ci_infCh[i] <- quantile(CDCh[i,], probs = 0.025)
  Ci_supCh[i] <- quantile(CDCh[i,], probs = 0.975)
}
CDstatsCh <- cbind(MeanCh, Ci_infCh, Ci_supCh)
rownames(CDstatsCh) <- rownames(CDCh)
colnames(CDstatsCh) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsCh, "cladeDensityStatsChi_100.csv")
write.csv(essimCh, "essimChi_100.csv")


CDDi <- matrix(ncol = length(trMa), nrow = nrow(overDi))
essimDi <- matrix(ncol = 2, nrow = length(trMa))
for(i in 1:length(trMa)) {
  claDenLogDi <- cladeDensity(mat = overDi, tree = trMa[[i]], log = T)
  CDDi[,i] <- claDenLogDi$CD[,2]
  essimDi[i,] <- claDenLogDi$essim[1:2]
}
rownames(CDDi) <- claDenLogDi$CD[,1]
colnames(essimDi) <- c("rho", "P")

MeanDi <- Ci_infDi <- Ci_supDi <- numeric()
for(i in 1:nrow(CDDi)) {
  MeanDi[i] <- mean(CDDi[i,])
  Ci_infDi[i] <- quantile(CDDi[i,], probs = 0.025)
  Ci_supDi[i] <- quantile(CDDi[i,], probs = 0.975)
}
CDstatsDi <- cbind(MeanDi, Ci_infDi, Ci_supDi)
rownames(CDstatsDi) <- rownames(CDDi)
colnames(CDstatsDi) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsDi, "cladeDensityStatsDip_100.csv")
write.csv(essimDi, "essimDip_100.csv")


CDPr <- matrix(ncol = length(trMa), nrow = nrow(overPr))
essimPr <- matrix(ncol = 2, nrow = length(trMa))
for(i in 1:length(trMa)) {
  claDenLogPr <- cladeDensity(mat = overPr, tree = trMa[[i]], log = T)
  CDPr[,i] <- claDenLogPr$CD[,2]
  essimPr[i,] <- claDenLogPr$essim[1:2]
}
rownames(CDPr) <- claDenLogPr$CD[,1]
colnames(essimPr) <- c("rho", "P")

MeanPr <- Ci_infPr <- Ci_supPr <- numeric()
for(i in 1:nrow(CDPr)) {
  MeanPr[i] <- mean(CDPr[i,])
  Ci_infPr[i] <- quantile(CDPr[i,], probs = 0.025)
  Ci_supPr[i] <- quantile(CDPr[i,], probs = 0.975)
}
CDstatsPr <- cbind(MeanPr, Ci_infPr, Ci_supPr)
rownames(CDstatsPr) <- rownames(CDPr)
colnames(CDstatsPr) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsPr, "cladeDensityStatsPri_100.csv")
write.csv(essimPr, "essimPri_100.csv")


CDSc <- matrix(ncol = length(trSq), nrow = nrow(overSc))
essimSc <- matrix(ncol = 2, nrow = length(trSq))
for(i in 1:length(trSq)) {
  claDenLogSc <- cladeDensity(mat = overSc, tree = trSq[[i]], log = T)
  CDSc[,i] <- claDenLogSc$CD[,2]
  essimSc[i,] <- claDenLogSc$essim[1:2]
}
rownames(CDSc) <- claDenLogSc$CD[,1]
colnames(essimSc) <- c("rho", "P")

MeanSc <- Ci_infSc <- Ci_supSc <- numeric()
for(i in 1:nrow(CDSc)) {
  MeanSc[i] <- mean(CDSc[i,])
  Ci_infSc[i] <- quantile(CDSc[i,], probs = 0.025)
  Ci_supSc[i] <- quantile(CDSc[i,], probs = 0.975)
}
CDstatsSc <- cbind(MeanSc, Ci_infSc, Ci_supSc)
rownames(CDstatsSc) <- rownames(CDSc)
colnames(CDstatsSc) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsSc, "cladeDensityStatsSci_100.csv")
write.csv(essimSc, "essimSci_100.csv")


CDGe <- matrix(ncol = length(trSq), nrow = nrow(overGe))
essimGe <- matrix(ncol = 2, nrow = length(trSq))
for(i in 1:length(trSq)) {
  claDenLogGe <- cladeDensity(mat = overGe, tree = trSq[[i]], log = T)
  CDGe[,i] <- claDenLogGe$CD[,2]
  essimGe[i,] <- claDenLogGe$essim[1:2]
}
rownames(CDGe) <- claDenLogGe$CD[,1]
colnames(essimGe) <- c("rho", "P")

MeanGe <- Ci_infGe <- Ci_supGe <- numeric()
for(i in 1:nrow(CDGe)) {
  MeanGe[i] <- mean(CDGe[i,])
  Ci_infGe[i] <- quantile(CDGe[i,], probs = 0.025)
  Ci_supGe[i] <- quantile(CDGe[i,], probs = 0.975)
}
CDstatsGe <- cbind(MeanGe, Ci_infGe, Ci_supGe)
rownames(CDstatsGe) <- rownames(CDGe)
colnames(CDstatsGe) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsGe, "cladeDensityStatsGek_100.csv")
write.csv(essimGe, "essimGek_100.csv")


CDIg <- matrix(ncol = length(trSq), nrow = nrow(overIg))
essimIg <- matrix(ncol = 2, nrow = length(trSq))
for(i in 1:length(trSq)) {
  claDenLogIg <- cladeDensity(mat = overIg, tree = trSq[[i]], log = T)
  CDIg[,i] <- claDenLogIg$CD[,2]
  essimIg[i,] <- claDenLogIg$essim[1:2]
}
rownames(CDIg) <- claDenLogIg$CD[,1]
colnames(essimIg) <- c("rho", "P")

MeanIg <- Ci_infIg <- Ci_supIg <- numeric()
for(i in 1:nrow(CDIg)) {
  MeanIg[i] <- mean(CDIg[i,])
  Ci_infIg[i] <- quantile(CDIg[i,], probs = 0.025)
  Ci_supIg[i] <- quantile(CDIg[i,], probs = 0.975)
}
CDstatsIg <- cbind(MeanIg, Ci_infIg, Ci_supIg)
rownames(CDstatsIg) <- rownames(CDIg)
colnames(CDstatsIg) <- c("Mean", "Ci_inf", "Ci_sup")
write.csv(CDstatsIg, "cladeDensityStatsIgu_100.csv")
write.csv(essimIg, "essimIgu_100.csv")


essimAg <- read.csv("essimAng_100.csv")
essimGe <- read.csv("essimGek_100.csv")
essimIg <- read.csv("essimIgu_100.csv")
essimSc <- read.csv("essimSci_100.csv")
essimCe <- read.csv("essimCet_100.csv")
essimCh <- read.csv("essimChi_100.csv")
essimDi <- read.csv("essimDip_100.csv")
essimPr <- read.csv("essimPri_100.csv")

pdf("Figure_4.pdf", width = 9, height = 6)
layout(matrix(1:8, ncol = 4, byrow = TRUE))
col1 <- rgb(70/255, 190/255, 173/255)
hist(essimAg$rho, breaks = 16, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Anguimorpha", xlab = expression(rho))
hist(essimGe$rho, breaks = 10, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Gekkota", xlab = expression(rho))
hist(essimIg$rho, breaks = 3, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Iguania", xlab = expression(rho))
hist(essimSc$rho, breaks = 4, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Scincoidea", xlab = expression(rho))
hist(essimCe$rho, breaks = 5, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Cetartiodactyla", xlab = expression(rho))
hist(essimCh$rho, breaks = 3, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Chiroptera", xlab = expression(rho))
hist(essimDi$rho, breaks = 5, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Diprotodontia", xlab = expression(rho))
hist(essimPr$rho, breaks = 4, xlim = c(-0.5,0.5), col = col1, border = NA, 
     main = "Primates", xlab = expression(rho))
dev.off()
