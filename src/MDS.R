#!/usr/bin/Rscript --slave

#Packages
#install.packages("ggfortify")
library(ggplot2)
library(ggfortify)
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
library(cluster)
#kohonen library
install.packages("kohonen")
library(kohonen)
library(MASS)

#Récupération des arguments
variable <- commandArgs(trailingOnly=TRUE)

#lecture du fichier csv
MyData <- read.csv(file="~/Bureau/projet_long/result/other/vecteur_db.csv", header=TRUE, sep=";",stringsAsFactors=FALSE)
#MyData <- read.csv(file=variable[1], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData[is.na(MyData)] <- 1 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData$ACC<-(MyData$ACC-mean(MyData$ACC))/sd(MyData$ACC)

MyData$polar<-(MyData$polar-mean(MyData$polar))/sd(MyData$polar)

MyData$hydrogenB<-(MyData$hydrogenB-mean(MyData$hydrogenB))/sd(MyData$hydrogenB)

MyData$hydrophobic<-(MyData$hydrophobic-mean(MyData$hydrophobic))/sd(MyData$hydrophobic)

MyData$vdw<-(MyData$vdw-mean(MyData$vdw))/sd(MyData$vdw)

MyData$ionic<-(MyData$ionic-mean(MyData$ionic))/sd(MyData$ionic)

########################"
MyData.Vector <- MyData[, 2:length(MyData)]
MyData.AA <- MyData$AA
par(mfrow=c(1,2))
autoplot(cmdscale(dist(MyData.Vector), eig = TRUE), label = TRUE, label.size = 3)
autoplot(isoMDS(dist(MyData.Vector)), shape = FALSE, label.colour = 'blue', label.size = 3)
################################################################# NMDS
sol <- metaMDS(MyData.Vector)
sol
plot(sol, type="t")

## Start from previous best solution
sol2 <- metaMDS(dune, previous.best = sol)


###############################################"

prox.mds <- isoMDS(dist(MyData.Vector),cmdscale(dist(MyData.Vector)))
X <- prox.mds$points[,1]
Y <- prox.mds$points[,2]
plot(prox.mds$points, pch=".")
text(X,Y,proxim[,2],pos=1,cex=0.8)


prox.sam <- sammon(dist(MyData.Vector),cmdscale(dist(MyData.Vector),k=2))
X <- prox.sam$points[,1]
Y <- prox.sam$points[,2]
plot(prox.sam$points,pch=".")
text(X,Y,proxim[,2],pos=1,cex=0.8)
##################################################""
require(vegan)

mod <- capscale(varespec ~ 1)
par(mfrow=c(1,2))
plot(mod)
plot(rda(varespec))

mod2 <- prcomp(varespec)
biplot(mod2)


library(ProjectionBasedClustering)

MDS(as.matrix(MyData.Vector),method='euclidean',OutputDimension=2,PlotIt=TRUE)
SammonsMapping(as.matrix(MyData.Vector),method= 'euclidean',OutputDimension=2,PlotIt=FALSE,Cls)


