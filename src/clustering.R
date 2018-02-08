#!/usr/bin/Rscript --slave

#Packages
#install.packages("ggfortify")
library(ggplot2)
library(ggfortify)
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
library(cluster)
#kohonen library
#install.packages("kohonen")
library(kohonen)
#install.packages(ProjectionBasedClustering)
library(ProjectionBasedClustering)

#Récupération des arguments
variable <- commandArgs(trailingOnly=TRUE)

#lecture du fichier csv
#MyData <- read.csv(file="~/Bureau/projet_long/result/other/vecteur_db.csv", header=TRUE, sep=";",stringsAsFactors=FALSE)
MyData <- read.csv(file=variable[1], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData[is.na(MyData)] <- 0 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData$ACC<-(MyData$ACC-mean(MyData$ACC))/sd(MyData$ACC)

MyData$polar<-(MyData$polar-mean(MyData$polar))/sd(MyData$polar)

MyData$hydrogenB<-(MyData$hydrogenB-mean(MyData$hydrogenB))/sd(MyData$hydrogenB)

MyData$hydrophobic<-(MyData$hydrophobic-mean(MyData$hydrophobic))/sd(MyData$hydrophobic)

MyData$vdw<-(MyData$vdw-mean(MyData$vdw))/sd(MyData$vdw)

MyData$ionic<-(MyData$ionic-mean(MyData$ionic))/sd(MyData$ionic)

##############################################################PCA
################################""ACP

MyData.Vector <- MyData[, 2:length(MyData)]
MyData.AA <- MyData$AA
MyData.Vector.pca <- prcomp(MyData.Vector)

a=variable[2]
pdf(paste(a,"cluster.pdf"))
plot(MyData.Vector.pca, type = "l",main="Variances en fonction des Composantes Principales")
summary(MyData.Vector.pca)




# Plot
# Colorer en fonction du cos2: qualité de représentation
fviz_pca_var(MyData.Vector.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)


################################MDS
MDS(as.matrix(MyData.Vector),method='euclidean',OutputDimension=2,PlotIt=TRUE)

################################SAMMON
#SammonsMapping(as.matrix(MyData.Vector),method= 'euclidean',OutputDimension=2,PlotIt=TRUE)

###############################################  #clusterisation
######################################## #kmeans
pourcentage <-function(n,v,y,z){
  polar=rownames(table(n, z$cluster))
  polar=lapply(polar, as.numeric)
  toto = table(n, z$cluster)
  test=matrix(nrow=1, ncol=dim(toto)[2])
  for (i in 1:dim(toto)[2])
  {
    mat= 0
    for(j in 1:dim(toto)[1]){
      coeff=toto[j,i]
      a=coeff*abs(as.numeric(polar[i]))
      mat=mat+a
      
    }
    test[i]=mat
  }

  sum_polar=sum(test[1,])

  for (i in 1:dim(test)[2]){
    test[,i]=(abs(as.numeric(test[,i]))*100)/abs(sum_polar)
  }
  colnames(test)=colnames(nb_obs)
  barplot(test,ylab=v,xlab="Cluster", main= paste(v,y),ylim=c(0,100))
  
  return(0)
}


    #DETERMINATION NB CLUSTER
nb_cluster=5
fviz_nbclust(MyData.Vector, kmeans, method = "wss")+
  geom_vline(xintercept = nb_cluster, linetype = 2)
    #CALCUL
set.seed(8953)
MyData2 <- MyData
MyData2$AA <- NULL
(kmeans.result <- kmeans(MyData2, nb_cluster))
    #plot
autoplot(kmeans.result, data = MyData, label = TRUE, label.size = nb_cluster, frame = TRUE,main="Clustering par Kmeans")

#calcul fréquences
#score=log2(fqobs/fqatt)
nb_obs=table(MyData$AA, kmeans.result$cluster)
fq_obs <- matrix(nrow=dim(nb_obs)[1], ncol=nb_cluster)
for (i in 1:dim(nb_obs)[1])
{
  for(j in 1:nb_cluster){
    a=nb_obs[i,j]/sum(nb_obs[,j])
    fq_obs[i,j]=a
  }
}
nb_att=table(MyData$AA)
fq_att=round(prop.table(nb_att), nb_cluster)#fréquence 

#matrice score par cluster
mat <- matrix(nrow=length(nb_att), ncol=nb_cluster)
for (i in 1:length(nb_att))
{
  for(j in 1:nb_cluster){
    a=log2(fq_obs[i,j]/fq_att[i])
    mat[i,j]=ifelse(a == -Inf,0,a)
  }
}
colnames(mat)=colnames(nb_obs)
rownames(mat)=rownames(nb_obs)
print("Matrice des scores par cluster:")
print(mat)


#par(mfrow=c(2,dim(mat)[2]/2))
for (i in 1:dim(mat)[2])
{
  barplot(mat[,i],main='Distribution des scores (Méthodes Kmeans)',xlab=paste("Cluster ",i),ylab="Score",ylim=c(4,-4))
}




######################################################  #kmedoid
    #Clustering with pam()


pam.res=pam(MyData.Vector, nb_cluster)
autoplot(pam.res, frame = TRUE, frame.type = 'norm',main="Clustering par Kmedoid(PAM)")
fviz_silhouette(silhouette(pam.res)) 
# Compute silhouette
sil <- silhouette(pam.res)[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]

#calcul fréquences
#log2(fqobs/fqatt)=score
nb_obs=table(MyData$AA, pam.res$cluster)
fq_obs <- matrix(nrow=dim(nb_obs)[1], ncol=nb_cluster)
for (i in 1:dim(nb_obs)[1])
{
  for(j in 1:nb_cluster){
    a=nb_obs[i,j]/sum(nb_obs[,j])
    fq_obs[i,j]=a
  }
}
nb_att=table(MyData$AA)
fq_att=round(prop.table(nb_att), nb_cluster)#fréquence 

#matrice score par cluster
mat <- matrix(nrow=length(nb_att), ncol=nb_cluster)
for (i in 1:length(nb_att))
{
  for(j in 1:nb_cluster){
    a=log2(fq_obs[i,j]/fq_att[i])
    mat[i,j]=ifelse(a == -Inf,0,a)
  }
}
colnames(mat)=colnames(nb_obs)
rownames(mat)=rownames(nb_obs)
print("Matrice des scores par cluster:")
print(mat)
#par(mfrow=c(2,dim(mat)[2]/2))
for (i in 1:dim(mat)[2])
{
  barplot(mat[,i],main='Distribution des scores (Méthode PAM)',xlab=paste("Cluster ",i),ylab="Score",ylim=c(4,-4))
}
########################################
pourcentage(MyData$ACC,"Accessibilité au solvant (%)","Kmeans",kmeans.result)
pourcentage(MyData$polar,"% d'interaction polar","Kmeans",kmeans.result)
pourcentage(MyData$ionic," % d'interaction ionic","Kmeans",kmeans.result)
pourcentage(MyData$hydrophobic," % d'interaction hydrophobic","Kmeans",kmeans.result)
pourcentage(MyData$vdw," % d'interaction van der Waals","Kmeans",kmeans.result)
pourcentage(MyData$hydrogenB," % de liaisons hydrogènes","Kmeans",kmeans.result)
#####################################
pourcentage(MyData$ACC,"Accessibilité au solvant (%)","Kmedoid",pam.res)
pourcentage(MyData$polar,"% d'interaction polar","Kmedoid",pam.res)
pourcentage(MyData$ionic," % d'interaction ionic","Kmedoid",pam.res)
pourcentage(MyData$hydrophobic," % d'interaction hydrophobic","Kmedoid",pam.res)
pourcentage(MyData$vdw," % d'interaction van der Waals","Kmedoid",pam.res)
pourcentage(MyData$hydrogenB," % de liaisons hydrogènes","Kmedoid",pam.res)


###########################################"kohonen"

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
plot(carte,type="count",palette.name=degrade.bleu,main=paste("Carte de Kohonen n=",n))

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



par(mfrow=c(1,1))

return(0)
}


test(4)

test(9)

test(16)

test(25)
#################################toroidal-----------------------------
thor<-function(n){
set.seed(100)
carte_tor <-som(as.matrix(MyData.Vector),grid=somgrid(xdim = n, ydim = n, topo =  "hexagonal",
neighbourhood.fct = "bubble", toroidal = TRUE))
#summary
print(summary(carte_tor))
#architecture of the grid
print(carte_tor$grid)

#count plot
plot(carte_tor,type="count",palette.name=degrade.bleu,main=paste("Carte de Kohonen toroidal n=",n))

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

return(0)
}

thor(4)

thor(9)

thor(16)

thor(25)
dev.off()