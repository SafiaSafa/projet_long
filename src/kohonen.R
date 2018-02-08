
#lecture du fichier csv



#kohonen library
install.packages("kohonen")
library(kohonen)

MyData <- read.csv(file="~/projet_long/result/other/vecteur_db.csv", header=TRUE, sep=";",stringsAsFactors=FALSE)
#MyData <- read.csv(file=variable[1], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData[is.na(MyData)] <- 0 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData$ACC<-(MyData$ACC-mean(MyData$ACC))/sd(MyData$ACC)

MyData$polar<-(MyData$polar-mean(MyData$polar))/sd(MyData$polar)

MyData$hydrogenB<-(MyData$hydrogenB-mean(MyData$hydrogenB))/sd(MyData$hydrogenB)

MyData$hydrophobic<-(MyData$hydrophobic-mean(MyData$hydrophobic))/sd(MyData$hydrophobic)

MyData$vdw<-(MyData$vdw-mean(MyData$vdw))/sd(MyData$vdw)

MyData$ionic<-(MyData$ionic-mean(MyData$ionic))/sd(MyData$ionic)

MyData.Vector <- MyData[, 2:length(MyData)]
MyData.AA <- MyData$AA
#####################################################
#SOM  Pour 4/9/16/25 groupe
degrade.bleu <-function(n){
    return(rgb(0,0.4,1,alpha=seq(0,1,1/n)))
}
coolBlueHotRed <-function(n, alpha = 1) {
    rainbow(n, end=4/6, alpha=alpha)[n:1]
}


thor<-function(n){
    set.seed(100)
    SelectVector<-MyData.Vector[, !apply(MyData.Vector == 0, 2, all)]
    carte_tor <-som(as.matrix(SelectVector),grid=somgrid(xdim = n, ydim = n, topo =  "hexagonal",
                                                         neighbourhood.fct = "bubble", toroidal = TRUE))
    #summary
    print(summary(carte_tor))
    #architecture of the grid
    print(carte_tor$grid)
    
    #count plot
    plot(carte_tor,type="count",palette.name=degrade.bleu,main=paste("Carte de Kohonen toroidal n=",n*n))
    
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
    #graphique pour chaque variable
    par(mfrow=c(7,4))
    for (j in 1:ncol(SelectVector)){
        plot(carte_tor,type="property",property=carte_tor$codes[[1]][,j],palette.name=coolBlueHotRed,main=colnames(SelectVector)[j],cex=0.5)
    }
    par(mfrow=c(1,1))
    return(0)
}

thor(2)

thor(3)

thor(4)

thor(5)
