# Projet Long: Méthode de génération de score de compatibilité 1D-3D par un algorithme d’apprentissage automatique

READ ME

*Auteur* : SAFA-TAHAR-HENNI Safia

*Date* : 09/02/2018


Ce programme permet de générer de score de compatibilité 1D-3D par un algorithme d’apprentissage automatique.


Il est nessecaire pour l'utiliser d'installer la suite d'outil PBxplore. Pour cela il faut suivre les instructions suivantes:
```{r, engine='bash',count_lines}
cd bin/
git clone https://github.com/pierrepo/PBxplore
cd PBxplore
pip2 install --user pbxplore
 pip2 install pytest
python2 setup.py test
```
Pour Arppegio vous aurez besoin des dépendances suivantes (https://bitbucket.org/harryjubb/arpeggio):

- Python (v2.7)

- Numpy

- BioPython (>= v1.60)

- OpenBabel (with Python bindings)


Packages R (v3.4.2) nécessaire à l'exécution de ce programme:

- rpart

- rpart.plot

- ggfortify

- ggplot2

- factoextra

- cluster

- kohonen

- ProjectionBasedClustering

Usage: 
```{r, engine='bash',count_lines}
python2 projet.py [-h] -i Directory_sequence_to_predict -db Directory_database -output OUTPUT_Directory
```


exemple: python2 src/scoring.py -i data/ -db data/DB -o Projet_long



Pour plus d'information:
	python2 scoring.py -h



Repertoire:

	- bin : dossier contenant la librairie PBxplore, arpeggio et l'executable Dssp
	
	- data : contient un exemple de fichier pdb pour la prediction et un base de donné contenant 10 pdb issue de HOMSTRAD
	
	- doc : pdf du rapport de projet long et de la présentation oral
	
	- Resultat : contient les résultats de l'utilisation du programe sur les données contenu dans data 
	
	- src : contient le programe principal scoring.py ainsi que les modules associés
