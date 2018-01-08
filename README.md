# projet long: Méthode de génération de score de compatibilité 1D-3D par un algorithme d’apprentissage automatique

READ ME

auteur: SAFA-TAHAR-HENNI Safia
date: 19/09/2017 


Ce programme permet de regrouper par les k-means de structures à l'aide d'un alphabet structural.

Il est nessecaire pour l'utiliser d'installer la suite d'outil PBxplore. Pour cela il faut suivre les instructions suivantes:
	git clone https://github.com/pierrepo/PBxplore
    	cd PBxplore
	pip2 install --user pbxplore
 	pip2 install pytest
	python2 setup.py test

Usage: 
	python projet.py [-h] --input_directory PATH_SEQUENCES --output OUTPUT_FILE [--ncluster NOMBRE_CLUSTER] [--step NOMBRE_PAS]

	exemple: python src/projet.py --input_directory data/DM_Calf-1_WT --output 2017_9_19_WT

Pour plus d'information:
	python projet.py -h

Repertoire:
	- bin : dossier contenant la librairie PBxplore
	- data : contient les resultats de simulations de dynamiques moléculaires venant du domaine Calf-1 (le dossier étant trop lourd pour Moodle les données sont diponible à l'adresse suivante : http://www.dsimb.inserm.fr/~debrevern/Calf-1Projekt/)
	- doc : documenation relative au script python 
	- Resultat : contient les résultats de l'utilisation du programe sur les données contenu dans data 
	- src : contient le programe principal projet.py ainsi que les modules associés
