
#lecture du fichier csv



#kohonen library
install.packages("kohonen")
library(kohonen)
#SOM  Pour 4/9/16/25 groupe
degrade.bleu <-function(n){
return(rgb(0,0.4,1,alpha=seq(0,1,1/n)))
}
coolBlueHotRed <-function(n, alpha = 1) {
rainbow(n, end=4/6, alpha=alpha)[n:1]
}


test<-function(n){
set.seed(100)
carte <-som(as.matrix(MyData.Vector),grid=somgrid(n,n,"hexagonal"))
#summary
print(summary(carte))
#architecture of the grid
print(carte$grid)


#count plot
plot(carte,type="count",palette.name=degrade.bleu)

#noeud d’appartenance des observations
print(carte$unit.classif)
#nombre d’observations affectés à chaque noeud
nb <-table(carte$unit.classif)
print(nb)
#check if there are empty nodes
print(length(nb))

#plot distance to neighbours
plot(carte,type="dist.neighbours")

#codebooks–profils des noeuds
plot(carte,type="codes",codeRendering = "segments")
#tableaudes codebooks pour les deux premiers noeuds
print(carte$codes[[1]][1:2,]) 

#graphique pour chaque variable
par(mfrow=c(6,4))
for (j in 1:ncol(MyData.Vector)){
plot(carte,type="property",property=carte$codes[[1]][,j],palette.name=coolBlueHotRed,main=colnames(MyData.Vector)[j],cex=0.5)
}
par(mfrow=c(1,1))

return(0)
}
test(2)
-----------------------------------toroidal-----------------------------

set.seed(100)
carte_tor <-som(as.matrix(MyData.Vector),grid=somgrid(xdim = 15, ydim = 10, topo =  "hexagonal",
neighbourhood.fct = "bubble", toroidal = TRUE))
#summary
print(summary(carte_tor))
#architecture of the grid
print(carte_tor$grid)



#jeu de couleurs pour les nœuds de la carte
degrade.bleu <-function(n){
return(rgb(0,0.4,1,alpha=seq(0,1,1/n)))
}
#count plot
plot(carte_tor,type="count",palette.name=degrade.bleu)

#noeud d’appartenance des observations
print(carte_tor$unit.classif)
#nombre d’observations affectés à chaque noeud
nb <-table(carte_tor$unit.classif)
print(nb)
#check if there are empty nodes
print(length(nb))

#plot distance to neighbours
plot(carte_tor,type="dist.neighbours")

#codebooks–profils des noeuds
plot(carte_tor,type="codes",codeRendering = "segments")
#tableaudes codebooks pour les deux premiers noeuds
print(carte_tor$codes[[1]][1:2,]) 


#jeu de couleurs pour le graphique
coolBlueHotRed <-function(n, alpha = 1) {
rainbow(n, end=4/6, alpha=alpha)[n:1]
}
#graphique pour chaque variable
par(mfrow=c(6,4))
for (j in 1:ncol(MyData.Vector)){
plot(carte_tor,type="property",property=carte$codes[[1]][,j],palette.name=coolBlueHotRed,main=colnames(MyData.Vector)[j],cex=0.5)
}
par(mfrow=c(1,1))

return(0)
}