#!/usr/bin/Rscript --slave


#Packages
# install.packages("rpart")
# install.packages("rpart.plot")
library("rpart")
library("rpart.plot")

#Récupération des arguments
variable <- commandArgs(trailingOnly=TRUE)

#lecture du fichier csv

MyData <- read.csv(file=variable[1], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData[is.na(MyData)] <- 0 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData$ACC<-(MyData$ACC-mean(MyData$ACC))/sd(MyData$ACC)

MyData$polar<-(MyData$polar-mean(MyData$polar))/sd(MyData$polar)

MyData$hydrogenB<-(MyData$hydrogenB-mean(MyData$hydrogenB))/sd(MyData$hydrogenB)

MyData$hydrophobic<-(MyData$hydrophobic-mean(MyData$hydrophobic))/sd(MyData$hydrophobic)

MyData$vdw<-(MyData$vdw-mean(MyData$vdw))/sd(MyData$vdw)

MyData$ionic<-(MyData$ionic-mean(MyData$ionic))/sd(MyData$ionic)

#Constitution groupe d'apprentissage et de test

MyData_train = MyData


#Construction de l'arbre de décision
target= AA ~H_ss+B_ss+E_ss+G_ss+I_ss+T_ss+S_ss+X._ss+o_ss+o_ss+o_ss+o_ss+o_ss+o_ss+o_ss+o_ss+ACC+a_PB+b_PB+c_PB+d_PB+e_PB+f_PB+g_PB+h_PB+i_PB+j_PB+k_PB+l_PB+m_PB+n_PB+o_PB+p_PB+polar+ionic+hydrophobic+vdw+hydrogenB  



a=variable[3]
pdf(paste(a,"predict.pdf"))
fq_att=table(MyData$AA)
vecteur=c(1/fq_att)
w = NULL#vecteur de poids depend des fréquences il doit mesurer 110

for (i in 1:dim(MyData_train)[1])
{
    w[i]=vecteur[which(names(vecteur)==MyData_train$AA[i])]
  
}

mytree <- rpart(target, data = MyData_train, method = "class", minsplit = 2, minbucket = 1, weights = w)
#minsplit est "le nombre minimum d'observations qui doivent exister dans un noeud pour qu'une tentative de division soit effectuée" et minbucket est "le nombre minimum d'observations dans un noeud terminal"
rpart.plot(mytree,main='Arbre de décision')

#Elagage
tree_ms3 = rpart(target, MyData_train, control = rpart.control(minsplit = 3))
tree_ms10 = rpart(target, MyData_train, control = rpart.control(minsplit = 10))

rpart.plot(tree_ms3, main = "Prunning: minsplit(nombre minimum d'observations dans un noeud) = 3")

rpart.plot(tree_ms10, main = "Prunning: minsplit(nombre minimum d'observations dans un noeud) =10")


#prediction sur le jeu test
#lecture du fichier csv

MyData_t<- read.csv(file=variable[2], header=TRUE, sep=";",stringsAsFactors=FALSE)

MyData_t[is.na(MyData_t)] <- 0 #remplace les NA (s'il y en as) par des 0

#normalisation par zscore
MyData_t$ACC<-(MyData_t$ACC-mean(MyData_t$ACC))/sd(MyData_t$ACC)

MyData_t$polar<-(MyData_t$polar-mean(MyData_t$polar))/sd(MyData_t$polar)

MyData_t$hydrogenB<-(MyData_t$hydrogenB-mean(MyData_t$hydrogenB))/sd(MyData_t$hydrogenB)

MyData_t$hydrophobic<-(MyData_t$hydrophobic-mean(MyData_t$hydrophobic))/sd(MyData_t$hydrophobic)

MyData_t$vdw<-(MyData_t$vdw-mean(MyData_t$vdw))/sd(MyData_t$vdw)

MyData_t$ionic<-(MyData_t$ionic-mean(MyData_t$ionic))/sd(MyData_t$ionic)

MyData_test = MyData_t
predictions = predict(mytree , MyData_test)

vec=NULL
for (i in 1:dim(predictions)[1])
{
  vec[i]=names(which(predictions[i,]==max(predictions[i,])))
}
print("Séquence prédites:")
print(vec)


dev.off()
