#!/usr/bin/Rscript --slave

#Packages
#install.packages("ggfortify")
library(ggplot2)
library(ggfortify)
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
library(cluster)

#Récupération des arguments
variable <- commandArgs(trailingOnly=TRUE)

#lecture du fichier csv
#MyData <- read.csv(file="../result/test/vecteur_db.csv", header=TRUE, sep=";",stringsAsFactors=FALSE)
MyData <- read.csv(file=variable[1], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData[is.na(MyData)] <- 0 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData$ACC<-(MyData$ACC-mean(MyData$ACC))/sd(MyData$ACC)

MyData$polar<-(MyData$polar-mean(MyData$polar))/sd(MyData$polar)

MyData$hydrogenB<-(MyData$hydrogenB-mean(MyData$hydrogenB))/sd(MyData$hydrogenB)

MyData$hydrophobic<-(MyData$hydrophobic-mean(MyData$hydrophobic))/sd(MyData$hydrophobic)

MyData$vdw<-(MyData$vdw-mean(MyData$vdw))/sd(MyData$vdw)

MyData$ionic<-(MyData$ionic-mean(MyData$ionic))/sd(MyData$ionic)

#PCA
  #center = TRUE,scale. = TRUE
MyData.Vector <- MyData[, 2:length(MyData)]
MyData.AA <- MyData$AA
MyData.Vector.pca <- prcomp(MyData.Vector)
#a="result/test/"
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

  #clusterisation
  #kmeans
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
fq_obs=round(prop.table(nb_obs), nb_cluster)#fréquence 
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
  barplot(mat[,i],main='Distribution des scores (Méthodes Kmeans)',xlab=paste("Cluster ",i),ylab="Score",ylim=c(0,-10))
}

  #kmedoid
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
fq_obs=round(prop.table(nb_obs), nb_cluster)#fréquence 
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
  barplot(mat[,i],main='Distribution des scores (Méthode PAM)',xlab=paste("Cluster ",i),ylab="Score",ylim=c(0,-10))
}

dev.off()