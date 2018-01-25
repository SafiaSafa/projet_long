#! /usr/bin/env python
# -*- coding: utf-8 -*-

#! /usr/bin/env python
# -*- coding: utf-8 -*-

import pbxplore as pbx
import random
import matplotlib.pyplot as plt
import os
import numpy as np
from string import ascii_lowercase
from collections import OrderedDict
from itertools import count,  izip
import math
import sys
import argparse
import time
import csv
import glob
from os.path import basename
import Bio.PDB

def definition_des_arguments():
    """
    Implemenation de la fonction 'help'.
    Verifie que les arguments --path_sequences et --output_file
    sont corrects.
    Renvoie la liste des arguments.
    """

    parser = argparse.ArgumentParser(
                                    description="Ce script permet un \
                                    creation fichier vecteur environnment \
                                    à partir d'un alphabet structural.")
    parser.add_argument('-i', action="store", type=str,
                        dest='path_to_predict', required=True, help="repertoire \
                        comptenant les donnees a predire")
    parser.add_argument('-db', action="store", type=str,
                        dest='path_db', required=True, help="repertoire \
                        comptenant les donnees qui serviront de base de données")
    parser.add_argument('-o',  dest='output_file', action="store",
                        type=str, required=True, help='Nom du repertoire \
                        où seront stockees les sorties')

    args = parser.parse_args()
    if not os.path.exists(args.path_db):
        sys.exit("Ce dossier "+args.path_db+" n'existe pas! \n")
    if not os.path.exists(args.path_to_predict):
        sys.exit("Ce dossier "+args.path_to_predict+" n'existe pas! \n")
    if os.path.exists(args.output_file):
        sys.exit("Ce dossier "+args.output_file+" existe deja! \n")
    return args

def liste_sequence(path_sequences):
    '''
    Parse le dossier contenant les fichiers pdb a analyses.
    Assigne a chaque PDB une sequence de blocs proteique.
    Retourne une liste de sequences de blocs proteique (sequences) ainsi
    que la liste des PDB (liste_PDB) presents dans le repertoire
    (path_sequences).
    '''
    sequences = []
    liste_PDB = []
    liste_chain = []
    for element in os.listdir(path_sequences):
        if element.endswith('.pdb'):
            liste_PDB.append(element)
            for chain_name, chain in pbx.chains_from_files([path_sequences+"/"+element]):
                liste_chain.append(chain_name)
                # Calcule angle phi et psi
                dihedrials = chain.get_phi_psi_angles()
                # assignation des blocks proteique
                pb_seq = pbx.assign(dihedrials)
                sequences.append(pb_seq)
    return sequences, liste_PDB, liste_chain


def vecteur_PB(liste_sequences_PB):
    """
        renvoie la conversion d'une séquence de proteines block en liste de vesteur binaire
    """
    aDict = dict(zip('abcdefghijklmnop',  range(0, 16)))
    liste_vecteur = []
    for j in range(len(liste_sequences_PB)):
        liste_v = []
        for i in range(len(liste_sequences_PB[j])):
            vecteur = np.zeros(16)
            vecteur = vecteur.astype(int)
            if liste_sequences_PB[j][i] != 'Z':
                vecteur[aDict[liste_sequences_PB[j][i]]]=1
                vecteur = ' '.join(str(x) for x in vecteur)
                liste_v.append((vecteur.split(),i+1))
            else:
                vecteur = ' '.join(str(x) for x in vecteur)
                liste_v.append((vecteur.split(),i+1))

        liste_vecteur.append(liste_v)

    return(liste_vecteur)


def Dssp(pdb,repertoire):
    aDict = dict(zip('HBEGITS-',  range(0, 8)))
    print("dico structure IIaire :",aDict)
    liste_pdb = []
    for element in glob.glob(pdb+"/*"):
        if element.endswith('.pdb'):
            os.system("./bin/dssp "+element+" -o "+repertoire+"/out.txt")
            os.system("sed -i -n '/#/,//p' " +repertoire+"/out.txt")
            f = open(repertoire+"/out.txt",'r')
            t = f.readlines()
            liste = []
            t=t[1:]
            for ele in t:
                aa=ele[13:14]
                ss=ele[16:17] 
                if (ss == ' ') or (ss == '-'): 
                    ss = '0'
                else:
                    ss = aDict[ss]
                acc=ele[35:38]
                liste.append((aa,ss,acc))
            os.system("rm "+repertoire+"/out.txt")
            liste_pdb.append(liste)
    return(liste_pdb)

def vecteur_DSSP(liste_sequences_dssp):
    liste_vecteur = []
    for j in range(len(liste_sequences_dssp)):
        for i in range(len(liste_sequences_dssp[j])):
            liste_v = []
            vecteur = np.zeros(16)
            vecteur = vecteur.astype(int)          
            vecteur[int(liste_sequences_dssp[j][i][1])]=1
            liste_v.append(liste_sequences_dssp[j][i][0])
            liste_v.append(tuple(vecteur))            
            liste_v.append(int(liste_sequences_dssp[j][i][-1]))
            # print(liste_v)
            liste_vecteur.append(liste_v)
    # print(liste_vecteur)
    return(liste_vecteur)


def parse_contact(file):
    f = open(file, "r")
    t = f.readlines()
    dico={}
    ligne=[]
    for ele in t:
        ligne.append(ele.split())


    for ele in ligne:
        if ele[0].split("/")[1] not in dico:
            dico[ele[0].split("/")[1]]={"vdw":0,"hydrogenB":0,"ionic":0,"hydrophobic":0,"polar":0}
        else:
            if ele[4]>0:
                dico[ele[0].split("/")[1]]['vdw']+=int(ele[4])
            if ele[6]>0:
                dico[ele[0].split("/")[1]]['hydrogenB']+=int(ele[6])
            if ele[10]>0:
                dico[ele[0].split("/")[1]]["ionic"]+=int(ele[10])
            if ele[13]>0:
                dico[ele[0].split("/")[1]]["hydrophobic"]+=int(ele[13])
            if ele[16]>0:
                dico[ele[0].split("/")[1]]["polar"]+=int(ele[16])
    return(dico)



def arpeggio(path_data):
    liste_arpeggio = []
    for element in glob.glob(path_data+"/*"):
        if element.endswith('.pdb'):
            os.system("python bin/arpeggio/arpeggio.py "+element)
            for ele in glob.glob(path_data+"/*"):
                if ele.endswith('.contacts'): 
                    liste_arpeggio.append(parse_contact(ele))
                if ele.endswith('.pdb'):
                    pass
                else:
                    os.system('rm '+ele)
    return(liste_arpeggio)


def write_csv(liste_vecteur, liste_dssp, liste_arpeggio, args, repertoir,file_name):
    '''
        
    '''
    adress_abs = os.getcwd()


    # create a csv file in the result directory
    fname = repertoir+file_name

    file = open(fname, "w")
    writer = csv.writer(file, delimiter=';')

    fieldnames = ['AA','H_ss','B_ss', 'E_ss', 'G_ss', 'I_ss','T_ss' , 'S_ss','-_ss','o_ss','o_ss','o_ss','o_ss','o_ss','o_ss','o_ss','o_ss','ACC','a_PB','b_PB','c_PB','d_PB','e_PB','f_PB','g_PB','h_PB','i_PB','j_PB','k_PB','l_PB','m_PB','n_PB','o_PB','p_PB','polar', 'ionic', 'hydrophobic', 'vdw', 'hydrogenB']
    writer.writerow(fieldnames)
    for i in range(len(liste_vecteur)):
        for j in range(len(liste_vecteur[i])):
            if str(liste_vecteur[i][j][-1]) in liste_arpeggio[i].keys():
                TUPLE=(tuple(liste_dssp[j][0])+liste_dssp[j][1]+(liste_dssp[j][-1],)+tuple(liste_vecteur[i][j][0])+tuple(liste_arpeggio[i][str(liste_vecteur[i][j][-1])].values()))
            else:
                TUPLE=(tuple(liste_dssp[j][0])+liste_dssp[j][1]+(liste_dssp[j][-1],)+ tuple(liste_vecteur[i][j][0])+tuple(['0','0', '0', '0', '0']))
            writer.writerow(TUPLE)
    pass

def init(path_data,args,repertoire,file_name):
    sequences, liste_PDB, liste_chain = liste_sequence(path_data)
    liste_vec = vecteur_PB(sequences)
    liste_dssp = Dssp(path_data,repertoire)
    liste_vec_dssp=vecteur_DSSP(liste_dssp)
    liste_arpeggio = arpeggio(path_data)
    write_csv(liste_vec, liste_vec_dssp,liste_arpeggio, args, repertoire,file_name)

    pass

if __name__ == '__main__':
    # Recuperation de moment du lancement du code
    debut = time.time()
    # Recuperation des arguments
    args = definition_des_arguments()
    adress_abs = os.getcwd()
    # Creation du repertoire qui va contenir les resultats
    repertoire = adress_abs+"/result/"+args.output_file
    os.system("mkdir "+repertoire)

    #Construction du fichier csv contenant les vecteurs associés à la base de données
    path_db = adress_abs+"/"+args.path_db
    init(path_db,args,repertoire,"vecteur_db.csv")
    # #Construction du fichier csv contenant les vecteurs associés à la séquences à prédire
    path_tp = adress_abs+"/"+args.path_to_predict
    init(path_tp,args,repertoire,"vecteur_to_predict.csv")

    # #script R
    #clustering
    os.system("./src/clustering.R "+repertoire+"/vecteur_db.csv "+repertoire+" >"+repertoire+"clustering.txt")
    # # #Prediction
    os.system("./src/predict.R "+repertoire+"/vecteur_db.csv "+repertoire+"/vecteur_to_predict.csv "+repertoire+" >"+repertoire+"predict.txt")
    # Recuperation de moment de la fin d'execution du code
    fin = time.time()
    print "Duree du programme : {:.2f} secondes".format(fin - debut)


# ./src/clustering.R result/test/vecteur_db.csv result/test/