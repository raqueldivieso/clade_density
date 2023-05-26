#### Calculate matrix of range overlaping between species of  Chiroptera 
### The same was done for Carnivora, Cetartiodactyla, Diprotodontia, Primates, and squamates (Anguimorpha, Gekkota, Iguania, Lacertoidea and Scincoidea.
### Geographical distributions polygons were obtained from available data on The IUCN Red List of Threatened Species database Version 2018-2 and .

setwd("C:/Users/raque/OneDrive - ufpr.br/Projetos_Paralelos/Range_overlap")

rm(list=ls())

library(rgdal)
library(rgeos)
library(ape)
library(stringr)
library(tidyr)
library(dplyr)
library(bigmemory)


#open the IUCN shapefile for all chirop species
sha<-readOGR("shapefiles/Chiroptera/data_0.shp")
invisible(spTransform(sha, CRS=CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))) #world cylindrical equal area reprojection #898spp
length(sha$BINOMIAL)

#removing the shape features classified as uncertain or introduced:
sha <- sha[sha$LEGEND != "Presence Uncertain",]
sha <- sha[sha$LEGEND != "Extant & Introduced (resident)",]
sha <- sha[sha$LEGEND != "Extant & Origin Uncertain (resident)",]
sha <- sha[sha$LEGEND != "Presence Uncertain & Origin Uncertain",]
sha <- sha[sha$LEGEND != "Extinct & Introduced",]
sha <- sha[sha$LEGEND != "Presence Uncertain & Vagrant",]
sha <- sha[sha$LEGEND != "Possibly Extant & Origin Uncertain (resident)",]
sha <- sha[sha$LEGEND != "Extinct & Origin Uncertain",]
sha <- sha[sha$LEGEND != "Probably Extant & Introduced (resident)",]
sha <- sha[sha$LEGEND != "Extant & Reintroduced (resident)",]
unique <- unique(sha$BINOMIAL)

#write.csv(unique, "Chiro_spp_to_phylo.csv") #write a list of species with IUCN geographic range data to build the phylogeny in http://vertlife.org/
#open the tree to use just spp with match in both spp list:
tr<-read.nexus("trees/chiroptera_tree.nex")
spp<-tr$tree_7361$tip.label
spp <- sort(spp)
spp<-str_replace_all(spp,"_"," ")

tt<-data.frame(crossing(var1 = 1:length(spp), var2 = 1:length(spp)))
list<-tt[ !duplicated(apply(tt, 1, sort), MARGIN = 2), ]
list<-as.data.frame(paste(list$var1, list$var2, sep="_"))
colnames(list)<-"comb"
list<-as.list(list)
list<-list$comb
tt<-NULL
tr<-NULL
unique<-NULL
list<-NULL

overQuant <- matrix(ncol = length(spp), nrow = length(spp))
rownames(overQuant) <- colnames(overQuant) <- spp

for (i in 1:length(liste)) {
  try({
    v<-data.frame(strsplit(liste[i], split='_', fixed=TRUE))
    x<-as.numeric(v[1,1])
    y<-as.numeric(v[2,1])
    a <- sha[sha$BINOMIAL == spp[x], ]
    b <- sha[sha$BINOMIAL == spp[y], ]
    if (all(is.na(over(a,b)))==FALSE){
      xx<-as.integer(round(gArea(gIntersection(a,b))*10000)) -> overQuant[x, y]
    } else {
      overQuant[x, y]<-NA}
     if( ((i %% 1000) == 0)==TRUE){
     write.csv(overQuant, "overlap_matrix/Chiroptera_overlap_matrix.csv")}
  },
  silent=TRUE)}


